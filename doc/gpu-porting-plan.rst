GPU porting plan (NVIDIA GH200, CSC Roihu)
==========================================

This page is a step-by-step plan for making GIMIC run on the GPU partition of
CSC's Roihu supercomputer (NVIDIA GH200 Grace Hopper nodes). It explains what
the current code does, why a line-by-line port would not pay off, which
restructuring gives the speed-up, and how to carry the work out in verifiable
steps. Nothing on this page has been implemented yet.

.. contents::
   :local:
   :depth: 2


Summary of the recommendation
-----------------------------

* **Do not rewrite the whole program.** Of the ~13,000 lines of GIMIC, only
  about 1,500 do numerical work (``bfeval.f90``, ``caos.f90``,
  ``jtensor.F90``, the grid loops in ``jfield.f90`` and ``integral.f90``).
  The rest is input parsing (``getkw``), grid construction, VTK/cube output,
  the Python wrapper, tools and tests. All of that is architecture-neutral
  and already validated by the test suite. Keep it.

* **Do rewrite the numerical core**, as a new *batched, tiled* kernel that
  evaluates the current-density tensor for thousands of grid points at once.
  This turns the per-point matrix-vector products (memory-bound, cannot use
  a GPU well) into matrix-matrix products (compute-bound, exactly what a
  GPU is built for) and lets the existing screening become sub-matrix
  selection, which gives near-linear scaling in molecule size.

* **Do the restructuring on the CPU first**, validate it against the existing
  tests, then offload the same kernel. The CPU version of the new kernel is
  itself a large speed-up for big molecules and becomes the reference for
  the GPU version.

* **Programming model:** Fortran + OpenACC directives + cuBLAS, compiled with
  ``nvfortran`` from the ``nvhpc/26.3`` module on Roihu. Reasons are given
  below; the alternatives (OpenMP offload, ``do concurrent``, CUDA Fortran,
  C++/CUDA, Python/CuPy) are discussed too.

* **Multi-GPU:** one MPI rank per GPU, grid tiles distributed round-robin.
  A GH200 node has 4 GPUs; the ``gpularge`` partition allows up to 10 nodes.

* **Effort:** roughly 3 to 4 months for one developer who knows Fortran and
  has some GPU experience, with a usable CPU speed-up after the first 4 to 6
  weeks.


What the current code does
--------------------------

The whole cost of GIMIC is in evaluating the 3x3 current-density tensor
:math:`\mathcal{J}(\mathbf r)` on grid points. For every point the code
(``jtensor.F90``, subroutine ``contract``) does the following, with
:math:`N` the number of (Cartesian) contracted basis functions:

1. Evaluate all basis functions :math:`\phi_\mu(\mathbf r)`, their gradients
   :math:`\partial_m \phi_\mu`, the GIAO derivatives
   :math:`\partial_{B_b} \phi_\mu = (\mathbf r \times \mathbf R_A)_b\,\phi_\mu`
   and the mixed second derivatives (``bfeval.f90``, ``caos.f90``). Cost
   :math:`O(N)` with distance-based screening per contraction.

2. Seven matrix-vector products with the :math:`N \times N` density
   matrices: :math:`D\phi`, :math:`P_b^{T}\phi` for :math:`b = x,y,z`, and
   :math:`D\,\partial_{B_b}\phi` for :math:`b = x,y,z` (``dgemv`` when built
   with BLAS, ``matmul`` otherwise). Cost :math:`14 N^2` flops and, more
   importantly, :math:`7 \cdot 8 N^2` bytes of memory traffic per point.

3. About twenty dot products of length :math:`N` to assemble the tensor.

The loops over grid points are OpenMP-parallel (``calc_jtensors``,
``integrate_current``, ``integrate_modulus``, ``integrate_acid``). MPI
support exists in ``parallel.F90`` but is incomplete: the integration
routines carry a ``TODO`` where the MPI reduction should be, and the Python
wrapper simply prepends a launcher string.

Why this structure is a poor fit for a GPU as it stands:

* A matrix-vector product reads the whole matrix once per point. For
  benzene/def2-TZVP (:math:`N = 252`, :math:`D` is 0.5 MB) the matrices sit
  in cache and the code is fine. For :math:`N = 5000` each point streams
  1.4 GB through the memory system, so a core spends its time waiting on
  memory rather than computing. A GPU offloaded point-by-point would be
  limited the same way, and the launch overhead of tiny kernels would
  dominate.

* The per-point routines walk pointer-linked derived types
  (``molecule_t -> atom_t -> basis_t -> contraction_t`` with pointer
  components, plus the ``bfeval_t`` scratch object). Fortran pointer
  components inside derived types are exactly what offload compilers handle
  worst (deep copies, no contiguity guarantees), so the data model has to
  be flattened before any directive can be applied.

* The screening (``filter_screened``) is decided per point and per
  contraction with a branch; on a GPU this must become a per-tile index
  list so that the linear algebra runs on dense sub-matrices.

Cost model. With :math:`n_p` grid points and :math:`N_{\text{act}}` basis
functions surviving screening near a point, the dominant cost is
:math:`14\,N_{\text{act}}^2\,n_p` flops:

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 30

   * - Case
     - :math:`N`
     - Points
     - Flops
   * - benzene test (3d, 30^3)
     - 252
     - 27,000
     - 2.4e10 (trivial, no GPU needed)
   * - nanographene, ~120 atoms, TZVP
     - ~2,500
     - 8e6 (200^3)
     - 7e14
   * - large system, ~400 atoms, TZVP
     - ~8,000
     - 8e6
     - 7e15 before screening

One H100 sustains roughly 30 to 50 TFLOP/s in FP64 ``dgemm``; the batched
formulation therefore turns the last case into minutes of GPU time. The same
work done as ``dgemv`` on a CPU is bounded by memory bandwidth, not flops,
and is one to two orders of magnitude slower per socket. This is why the
reformulation, not the offload itself, is the heart of the plan.


The target: Roihu GPU partition
-------------------------------

Facts taken from the CSC documentation (Roihu system page, compiling page
and partition page; verify before use, the machine is new and the
documentation changes):

* 132 GPU nodes, each with 4 NVIDIA GH200 superchips: one H100 GPU with
  96 GiB HBM3 and one 72-core Grace (Arm Neoverse V2) CPU with 120 GiB
  LPDDR5 per superchip. CPU and GPU are cache-coherent over NVLink-C2C, so
  the GPU can address host memory directly (this is the "unified memory"
  that GH200 is known for).

* The GPU login node ``roihu-gpu.csc.fi`` is Arm. The CPU partition (AMD
  Turin, x86) and its login node are a different architecture: **binaries
  must be built on the GPU login node** for the GPU nodes. This also means
  the Python side of GIMIC (the ``gimic`` wrapper, Cython modules) needs an
  aarch64 Python environment there.

* Toolchains: default ``gcc/14.3.0 cuda/12.9.1 openmpi/5.0.10
  openblas/0.3.30``; alternatively ``gcc/15.2.0 cuda/13.1.1``; and
  ``nvhpc/26.3`` which bundles CUDA, an MPI and BLAS. GPU compute capability
  is 9.0 (``-gpu=cc90`` for nvhpc, ``-gencode arch=compute_90a,code=sm_90a``
  for nvcc). nvhpc supports OpenACC (``-acc=gpu``), OpenMP offload
  (``-mp=gpu``) and, for C++, ``-stdpar=gpu``.

* Partitions: ``gputest`` (15 min, 2 nodes), ``gpumedium`` (36 h, 1 node,
  1 to 4 GPUs), ``gpularge`` (36 h, up to 10 nodes, needs a scaling test),
  ``gpuinteractive`` with MIG slices (1/7 of a GPU, 12 GiB) planned. Each
  reserved GPU comes with 72 CPU cores and about 212 GiB of memory.

Consequences for the design:

* FP64 throughput on Hopper is high (tensor-core ``dgemm``), so there is no
  need to consider single precision. Keep FP64 everywhere: the diamagnetic
  and paramagnetic contributions cancel to a large degree and single
  precision would corrupt the small remainder.

* 96 GiB of HBM holds the four :math:`N \times N` FP64 matrices
  (:math:`D, P_x, P_y, P_z`) up to :math:`N \approx 50{,}000` (open-shell:
  :math:`N \approx 35{,}000`). Larger cases can leave the full matrices in
  the 120 GiB Grace memory and gather only the screened sub-blocks; the C2C
  link makes that practical.

* MIG slices make small jobs (integration planes, the "squares profile"
  workflow with hundreds of small runs) cheap once ``gpuinteractive`` is
  configured. The port should therefore work well on a fraction of a GPU
  too, which argues for modest per-tile memory use.


Design decisions and why
------------------------

Batched formulation of the tensor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a tile of :math:`n_p` points, build the matrices (points fastest in
memory)

.. math::

   \Phi_{p\mu} = \phi_\mu(\mathbf r_p), \qquad
   G^{(m)}_{p\mu} = \partial_m \phi_\mu(\mathbf r_p), \quad m = x,y,z .

Then two GEMM calls give everything the tensor needs:

.. math::

   [\,Y^{(0)} \;|\; Y^{(x)} \;|\; Y^{(y)} \;|\; Y^{(z)}\,] =
   [\,\Phi \;|\; G^{(x)} \;|\; G^{(y)} \;|\; G^{(z)}\,]\; D ,
   \qquad
   Z^{(b)} = \Phi\, P_b, \quad b = x,y,z .

(With the point-fastest layout, ``Y = Phi * D`` computes :math:`D^{T}\phi_p`
for every point, which is what the present code computes; :math:`D` is
symmetric. ``Z = Phi * P_b`` computes :math:`P_b^{T}\phi_p`, matching the
``dgemv('t', ...)`` in ``contract``. :math:`P_b` is antisymmetric, so this
transpose determines the sign of the paramagnetic term. Keep it.)

The tensor at point :math:`p` is then a per-point reduction over
:math:`\mu` that needs no further matrix operations. With
:math:`\mathbf d^{(a)} = \mathbf r_p \times \mathbf R_a` (the ``dbop`` of
``bfeval.f90``), :math:`a(\mu)` the atom of function :math:`\mu`, and
:math:`\varepsilon` the Levi-Civita symbol:

.. math::

   \rho_{\text{dia}} = \sum_\mu \Phi_{p\mu} Y^{(0)}_{p\mu}

   \text{ppd}_{mb} = \sum_\mu Z^{(b)}_{p\mu} G^{(m)}_{p\mu}

   \text{prsp1}_{mb} = -\sum_\mu d^{(a(\mu))}_b\, \Phi_{p\mu}\, Y^{(m)}_{p\mu}

   \text{prsp2}_{mb} = \sum_\mu Y^{(0)}_{p\mu}
        \big[\, d^{(a(\mu))}_b\, G^{(m)}_{p\mu}
             + \varepsilon_{bmd} R_{a(\mu),d}\, \Phi_{p\mu} \big]

   \mathcal J_{mb} = \zeta\,(\text{ppd} + \text{prsp1} + \text{prsp2})_{mb}
        + \text{diamagnetic off-diagonal terms } \pm \tfrac12 r_b\,\rho_{\text{dia}}

This is algebraically identical to ``contract`` (``prsp1`` uses the symmetry
of :math:`D` to move it from the :math:`\partial_B\phi` side to the
:math:`\partial_r\phi` side; ``prsp2`` expands the ``d2fdrdb`` array in
place). The explicit ``dbvec`` and ``d2fvec`` arrays disappear. Only the
atom index and coordinates per basis function are needed, which are two
small arrays.

Spin cases (``total``, ``alpha``, ``beta``, ``spindens``) only change which
:math:`D, P_b` are used; :math:`\Phi, G` are shared and should be evaluated
once per tile for all requested spin cases.

The ``spherical`` option (``cao2sao``) is best handled once at start-up by
transforming the density matrices to the Cartesian basis,
:math:`D_{\text{cart}} = C^{T} D_{\text{sph}} C`, so the kernel is always
Cartesian and pays nothing at run time.

Spatial tiles and screening
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Group grid points into spatially compact tiles (for the structured grids:
blocks of, say, 8x8x8 to 16x16x16 points; for ``Grid(file)`` grids: sort
the points along a Morton/Z-order curve and cut into chunks). For each tile
compute its bounding sphere (centre :math:`\mathbf c`, radius
:math:`\varrho`) and the *active set* :math:`A` of contractions with
:math:`|\mathbf c - \mathbf R_a| - \varrho \le \text{thrs}_i`, using the
same per-contraction radii that ``setup_screening`` already computes. Then

* :math:`\Phi, G` are only evaluated and stored for :math:`\mu \in A`
  (:math:`n_p \times N_{\text{act}}`),
* the GEMMs use the gathered sub-matrices :math:`D_{AA}` and
  :math:`P_{b,AA}` (gathered on the GPU, where the full matrices live),
* the reductions run over :math:`A` only.

For a large molecule :math:`N_{\text{act}}` saturates at the number of
functions within the screening radius of the tile, so the cost per point
becomes independent of molecule size. The gather costs
:math:`O(N_{\text{act}}^2)` per tile and is amortised over the
:math:`n_p` points of the tile, so tiles should hold a few thousand points.
Screening off (``Advanced.screening=off``) is simply :math:`A =` all.

Programming model
^^^^^^^^^^^^^^^^^

Recommended: **Fortran + OpenACC + cuBLAS with nvfortran**.

* The two custom kernels (basis-function evaluation, per-point reduction)
  are simple, embarrassingly parallel loops over ``(point, contraction)``
  and ``(point)``; they are not where the flops are. Directives express
  them adequately and keep the code readable to the current developers.
* The flops are in ``dgemm``, which is a cuBLAS call. nvhpc ships a Fortran
  ``cublas`` module; ``!$acc host_data use_device`` passes device pointers.
* OpenACC in nvfortran is the most mature directive path for Fortran on
  NVIDIA hardware, has explicit data-region control (which the tile pipeline
  needs), async queues for overlapping tiles, and ``acc_set_device_num`` for
  one-rank-per-GPU MPI.
* The same source compiles without the directives with gfortran (they are
  comments), so there is one code path for CPU and GPU, and the test suite
  runs on both.

Alternatives, and why they are not the first choice:

* **OpenMP target offload.** Equivalent in expressiveness; nvhpc supports it
  (``-mp=gpu``). Keeping both directive sets doubles maintenance; choose one.
  If the group prefers a single standard for CPU threading and GPU offload,
  OpenMP is acceptable, and because all directives will sit in one module
  (see Step 4) switching later is a bounded job.
* **Standard-parallel Fortran** (``do concurrent`` with ``-stdpar=gpu``).
  Zero directives and very attractive on GH200 because it relies on
  unified memory. It is a reasonable choice for the two custom kernels, but
  cuBLAS interop, explicit residency of the density matrices and async
  overlap are clumsier. Worth revisiting once nvhpc's stdpar Fortran
  matures; not the first choice today.
* **CUDA Fortran.** Full control, but nvfortran-only and more code. Not
  needed since the hot spot is cuBLAS.
* **Full rewrite in C++/CUDA.** Gives nothing over the above for a
  GEMM-dominated kernel, and would force rewriting the parser, grid, output
  and Python glue (about 90 % of the code) that gain nothing from being in
  C++.
* **Full rewrite in Python + CuPy.** The batched formulation is ~200 lines
  of NumPy/CuPy and would reach nearly the same cuBLAS-bound performance.
  This is the right choice only if the group wants to leave Fortran
  entirely; then the MOL/XDENS readers, grids and VTK writers must be
  re-implemented (``pygimic`` is a partial start), and deployment becomes a
  Python-environment problem on Arm. It is, however, an excellent tool for
  prototyping the formulas (see Step 2).

Data residency and the GH200 unified memory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Keep :math:`D, P_b` (and the flattened basis data) resident on the GPU for
the whole run with ``!$acc enter data copyin``. Per tile, the only host to
device traffic is the tile's point coordinates and active-set index list,
and device to host is :math:`9\,n_p` doubles of tensors. Explicit data
management is preferred over ``-gpu=mem:unified`` so the code also runs on
non-coherent GPUs (Puhti/Mahti/LUMI-style V100/A100/MI250 with OpenACC via
other compilers) and so that performance does not depend on page migration
heuristics. Unified memory remains a convenient debugging mode
(``-gpu=mem:unified`` compiles the same source).

Parallelism across GPUs and nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One MPI rank per GPU (``srun --ntasks-per-node=4 --gpus-per-node=4``),
each rank binding to ``mod(rank, ndevices)`` and using its 72 Grace cores
with OpenMP for tile preparation and output. Tiles are distributed
round-robin (tile :math:`t` to rank :math:`t \,\mathrm{mod}\, n`) which balances the
load under screening better than the contiguous split in ``schedule``.
Results are ``MPI_Gatherv``'d (tensors) or ``MPI_Reduce``'d (integrals) to
rank 0, which alone writes output, as now.

For the profile/"squares" workflows that launch hundreds of small
independent GIMIC runs, no MPI is needed: keep launching independent jobs
and give each one a GPU or a MIG slice.


Step-by-step guide
------------------

Each step ends in a state where the test suite passes. Steps 1 to 3 are
pure CPU work and give a large speed-up on their own for big molecules;
they are also the prerequisite for any GPU code. Effort estimates assume
one developer.

Step 0: Baseline and toolchain on Roihu (1 week)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. On ``roihu-gpu.csc.fi``: ``module load nvhpc/26.3``, create an aarch64
   Python environment with the ``requirements.txt`` packages (cython, numpy,
   runtest 2.3.2). Build GIMIC CPU-only with ``nvfortran``::

     ./setup --fc=nvfortran --cc=nvc --cxx=nvc++ --omp --blas \
             --extra-fc-flags="-O2 -Mpreprocess" build-nvhpc

   Fix whatever does not compile. Known risks: the GNU extensions
   ``etime``, ``dtime``, ``fdate``, ``hostnm`` (``timer.f90``,
   ``jfield.f90``, ``gimic.F90``, ``parallel.F90``); nvfortran supports
   most of them but check. Replace ``etime``/``dtime`` with
   ``system_clock`` while at it. ``.f90`` files that use preprocessor
   macros must be renamed ``.F90`` or compiled with ``-Mpreprocess``.
   CMake identifies nvfortran as compiler id ``NVHPC``; make sure
   ``cmake/downloaded/autocmake_fc.cmake`` and the OpenMP/BLAS modules do
   not assume GNU or Intel.

2. Run ``ctest`` with the nvfortran build and with the gcc/openblas build;
   both must pass. Record run times.

3. Create two or three *benchmark* inputs that represent the real target
   workloads (a 100 to 400 atom system with TZVP quality basis, a 3D
   ``cdens`` grid with :math:`10^6` to :math:`10^7` points, a
   ``Grid(bond)`` integration). The benzene tests are too small to measure
   anything. Time the baseline with 1, 72 and 288 OpenMP threads.

4. Add a lightweight timing report to ``gimic.bin`` (time in basis
   evaluation, contraction, output) so later steps can be measured without
   a profiler.

Why first: nothing else can be validated without a working nvfortran build
and realistic benchmarks, and compiler-portability problems are cheapest to
find in unmodified code.

Step 1: Flatten the data model (1 to 2 weeks)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add a module (``basis_flat.f90``) that, after ``new_basis`` and
``new_dens``, builds plain contiguous arrays from the pointer-linked types:

* per contraction (shell): ``l``, ``npf``, primitive offset, atom index,
  first Cartesian function index, ``ncomp``, screening radius;
* per primitive: exponent, normalised coefficient (``ncc``);
* per atom: coordinates;
* per Cartesian function: atom index, exponents :math:`(n_x, n_y, n_z)`
  (from ``gtodefs``), shell index;
* density matrices as ``real(8), allocatable :: D(N,N), P(N,N,3)`` per spin,
  already in the Cartesian basis if ``spherical`` is on (apply the
  ``cao2sao`` operator once), already reordered for Turbomole input.

Keep the old types for input parsing and output. Do not touch the kernels
yet.

Why: allocatable, contiguous arrays with no pointer components are the only
data that OpenACC/OpenMP offload can move and index efficiently, and this
step is independent of any GPU decision.

Step 2: Batched CPU kernel (2 to 3 weeks)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Write ``jtensor_batch.f90`` with one public entry point::

   subroutine jtensor_tile(kern, npts, r, spincase, T)
     ! r(3,npts) in, T(9,npts) out; kern holds the flat basis and densities

implemented as: (a) evaluate ``Phi(npts,N)``, ``G(npts,N,3)`` with a loop
over points and shells (each shell: one exponential sum over primitives,
then the Cartesian polynomial factors, mirroring ``cgto``/``dcgto``;
compute both value and gradient from the same exponentials);
(b) two ``dgemm`` calls (OpenBLAS/MKL/NVPL); (c) the per-point reduction
above.

Validation: a unit test that compares ``jtensor_tile`` against the old
``ctensor`` point by point on the benzene MOL/XDENS for random points,
including GIAO on/off, diamag/paramag on/off, open-shell, and spherical.
Agreement should be at the :math:`10^{-12}` relative level (different
summation order only). Then switch ``calc_jtensors`` and the three
``integrate_*`` routines to call the tile kernel on chunks of points and
run ``ctest``; loosen ``runtest`` tolerances only if a difference is
understood.

Optional but recommended before writing Fortran: prototype the formulas in
NumPy (about 200 lines, plus a small reader for MOL/XDENS) and check them
against ``jvec.vti`` from the benzene ``3d`` test. It costs a day and
removes the risk of chasing a sign error in Fortran later.

Why before GPU: this step contains all the algorithmic risk (index
bookkeeping, transposes, signs) and none of the GPU risk. It is also where
most of the CPU speed-up for large :math:`N` comes from, because ``dgemm``
reads the density matrix once per tile instead of once per point.

Step 3: Tiling and screening (1 to 2 weeks)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Tile builder: for structured grids, blocks in ``(i,j,k)``; for file
   grids, Morton sort then chunk. Each tile: point list, centre, radius.
2. Active set per tile from the screening radii; map to Cartesian function
   indices; gather ``D_AA``, ``P_AA(:,:,3)``.
3. ``jtensor_tile`` works on ``A`` only.
4. Choose the tile size adaptively so that ``Phi``/``G``/``Y``/``Z``
   (11 arrays of ``npts x N_act`` doubles) fit a memory budget, e.g. 1 to
   2 GB on the GPU later, a few hundred MB per thread on the CPU.
5. Test: results must be identical to Step 2 with screening off, and within
   the screening threshold with it on; the existing ``keyword-radius`` and
   ``giao-test`` cases cover the relevant paths.

Why: this converts the per-point branchy screening into dense sub-matrix
work and delivers the linear scaling in molecule size, which is the second
largest win after batching, and it defines the unit of work that the GPU
pipeline and the MPI distribution will use.

Step 4: GPU offload of the tile kernel (3 to 4 weeks)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All directives live in ``jtensor_batch.f90`` and a small ``gpu_env.f90``
(device selection, cuBLAS handle, memory budget query). Nothing outside
these files knows about the GPU.

1. Residency: after Step 1, ``!$acc enter data copyin`` the flat basis
   arrays and ``D``, ``P`` for every spin. Allocate the per-tile work arrays
   once at the maximum tile size (``!$acc enter data create``) and reuse.
2. Basis kernel: ``!$acc parallel loop gang vector collapse(2)`` over
   ``(shell, point)``, each iteration writing ``ncomp`` values of ``Phi``
   and ``G``; layout ``(point, function)`` so writes are coalesced across
   points. Exponentials are computed once per shell and point.
3. Gather kernel: ``D_AA(i,j) = D(A(i),A(j))`` on the device.
4. GEMMs: ``use cublas``; inside ``!$acc host_data use_device(Phi,G,DAA,Y)``
   call ``cublasDgemm`` (one call for ``Y``, one per ``b`` for ``Z``, or
   stack the three ``P_b`` into one ``N x 3N`` matrix for a single call).
   cuBLAS uses the FP64 tensor cores automatically.
5. Reduction kernel: ``!$acc parallel loop`` over points, inner sequential
   loop over ``A`` with the 12 scalar accumulators (9 tensor components
   plus the diamagnetic scalar and helpers); reads are coalesced across
   points thanks to the layout. Write ``T(9, npts)``.
6. Copy ``T`` back (``!$acc update self``) into the caller's array.
7. Overlap: two tile buffers on two async queues; while the GPU processes
   tile ``t``, the host builds the active set of ``t+1``. On GH200 the
   copies are small; the overlap mainly hides the tile bookkeeping.
8. Build flags: ``-acc=gpu -gpu=cc90 -cudalib=cublas -Minfo=accel``.
   Add CMake options ``ENABLE_OPENACC`` and ``ENABLE_CUBLAS`` that add these
   flags and link; without them the same source builds a CPU-only binary.
9. Validate against Step 3 on the benchmarks and the test suite. Profile
   with ``nsys`` to confirm the GEMMs dominate; use ``ncu`` only if the
   basis or reduction kernels take more than about 20 % of the time.

Expected outcome: a single H100 should process on the order of :math:`10^6`
points per minute for :math:`N_{\text{act}} \approx 2{,}000` to
:math:`3{,}000`, limited by ``dgemm``. The exact number depends on the
screening and is what Step 8 measures.

Step 5: Drivers, spin cases and integration on the tile API (1 week)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``calc_jtensors`` becomes a loop over tiles that writes into
  ``this%tens``; ``compute_jvectors``, ``jvector_plots`` and the VTK
  writers are unchanged.
* ``integrate_current``, ``integrate_modulus`` and ``integrate_acid`` are
  reduced to one routine: obtain ``T`` per tile, then apply the quadrature
  weights and the positive/negative, modulus and ACID bookkeeping on the
  host with OpenMP (the per-point post-processing is cheap; keeping it on
  the host keeps the GPU code small). The three-level weight structure of
  the Gauss grids maps onto the tile's ``(i,j,k)`` indices.
* Open-shell runs request ``total``, ``alpha``, ``beta`` and ``spindens``;
  compute all four from one evaluation of ``Phi``/``G`` per tile instead
  of four passes over the grid.
* Replace ``jfield_eta`` (which times 100 calls of the old kernel) by timing
  the first tile.
* Remove the old point-wise kernel once everything passes, or keep it
  behind a ``--reference`` flag for a release cycle. Keeping two kernels
  forever is not recommended.

Step 6: MPI across GPUs and nodes (1 to 2 weeks)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. ``init_mpi``: after ``MPI_Init``, ``acc_set_device_num(mod(rank,
   acc_get_num_devices(acc_device_nvidia)), acc_device_nvidia)``.
2. Distribute tiles round-robin; ``MPI_Gatherv`` tensors (``cdens``) and
   ``MPI_Reduce`` the integral accumulators (fix the ``TODO`` in
   ``integral.f90``). Only rank 0 writes files (already the case).
3. Use the MPI bundled with nvhpc (``mpif90`` wrapper); it is CUDA-aware but
   that is not required, as all reductions are on host arrays.
4. Slurm script for Roihu (check the exact ``--gpus`` syntax in the CSC
   documentation)::

     #!/bin/bash
     #SBATCH --account=project_XXXXXXX
     #SBATCH --partition=gpumedium
     #SBATCH --nodes=1
     #SBATCH --ntasks-per-node=4
     #SBATCH --cpus-per-task=72
     #SBATCH --gpus-per-node=4
     #SBATCH --time=02:00:00
     module load nvhpc/26.3
     export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
     srun gimic.bin gimic.inp > gimic.out

   The Python wrapper already accepts a launcher string; either pass
   ``srun`` through it or call ``gimic.bin`` directly as above.
5. Scaling test on ``gputest``/``gpumedium`` (1, 2, 4 GPUs) and then the
   multi-node run required to be allowed on ``gpularge``.

Step 7: Build, packaging and documentation (1 week)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``setup``/CMake: options ``--openacc``, ``--cublas``, ``--gpu-arch=cc90``;
  compiler-id handling for ``NVHPC``; a CPU-only build must still work with
  gfortran and OpenBLAS, and ``ctest`` must pass in both.
* An Apptainer definition based on the ``nvcr.io/nvidia/nvhpc`` image
  (``container/`` already has a CPU recipe) so the Arm build is
  reproducible; run it with ``--nv`` on the GPU nodes.
* CI: GitHub Actions cannot test the GPU path; add an nvfortran CPU-only
  job (the nvhpc container is public) so the directive code at least
  compiles in CI, plus the existing gfortran job.
* Documentation: installation on Roihu, the new options, tile-size and
  memory settings, expected accuracy, and this page updated to "as built".

Step 8: Validation and acceptance (1 week, overlaps with the above)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Full ``ctest`` on: gfortran/x86 CPU, nvfortran/Arm CPU-only, nvfortran
  GPU, and 4-rank MPI+GPU. Same reference files for all four.
* Benchmarks from Step 0: report time and energy per :math:`10^6` points
  for CPU-node (batched kernel) and 1/4 GPUs, and the speed-up over the
  current code. Expect the GPU to matter above roughly :math:`N = 1{,}000`
  and :math:`10^5` points; below that the CPU batched kernel is enough and
  should be what users run on the CPU partition.
* Numerical acceptance: integrated currents agree with the old code to
  better than :math:`10^{-6}` nA/T; vector fields agree pointwise to
  :math:`10^{-10}` relative.


Optional later optimisations
----------------------------

In decreasing order of expected value:

* Skip GEMM columns for functions whose values in a tile are all below
  threshold (a second, cheaper screening after evaluation).
* Exploit the symmetry of :math:`D` and antisymmetry of :math:`P_b` with
  ``dsymm``-type calls (at most a factor 2 on the GEMM).
* Several magnetic-field directions from one evaluation of :math:`\Phi`
  (the tensor is field-independent; only the final contraction with
  :math:`\mathbf B` changes). The current ``magnet``/``magnet_axis``
  design already separates these.
* Evaluate primitives in FP32 for the *screening estimate* only; never for
  the actual values.
* Sparse ``P_b`` for very large molecules with sparse perturbed densities.


Risks and how the plan contains them
------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Risk
     - Mitigation
   * - nvfortran rejects GNU extensions or old constructs
     - Step 0 builds unmodified code first; fixes are small and benefit
       every compiler.
   * - Sign or transpose error in the batched formulas
     - Step 2 unit test against the old kernel point by point, before any
       GPU work; optional NumPy prototype.
   * - Screening in tiles changes results
     - Step 3 compares screening on/off and to the old per-point screening;
       radii are the conservative ones already in ``basis.f90``.
   * - GPU memory for very large :math:`N`
     - Full matrices may stay in Grace memory; tiles use gathered
       sub-blocks; tile size adapts to a memory budget.
   * - Small jobs do not fill a GPU
     - Tile size and MIG slices; the batched CPU kernel remains the right
       tool for benzene-sized runs and for the CPU partition.
   * - Toolchain drift on a new machine
     - Keep the CPU-only build path first-class; container recipe pins
       nvhpc.
   * - Only one developer knows the GPU code
     - All GPU code in two files behind one Fortran interface; the rest of
       GIMIC is untouched.


Appendix: pseudo-code of the tile kernel
----------------------------------------

::

   for tile in tiles(grid):                          # host, OpenMP
       A  = active_functions(tile)                    # from screening radii
       upload tile.points, A                          # small
       gather  DAA = D(A,A), PAA(:,:,b) = P(A,A,b)    # device kernel
       Phi, G  = eval_basis(tile.points, A)           # device kernel
       Y       = [Phi G1 G2 G3] * DAA                 # cublasDgemm
       Z_b     = Phi * PAA(:,:,b),  b = 1..3          # cublasDgemm
       T       = reduce(Phi, G, Y, Z, atom(A), R)     # device kernel
       download T(9, npts)                            # 72 bytes/point
       post-process on host: J = T.B, weights, sums, ACID, VTK buffers

Memory per tile with :math:`n_p` points and :math:`N_{\text{act}}` active
functions: :math:`11 \cdot 8 \cdot n_p N_{\text{act}}` bytes for the work
arrays plus :math:`4 \cdot 8 \cdot N_{\text{act}}^2` for the sub-matrices.
With :math:`n_p = 4096` and :math:`N_{\text{act}} = 3000` this is about
1.1 GB plus 0.3 GB, comfortable on a 96 GiB H100 and workable on a 12 GiB
MIG slice with smaller tiles.


Sources
-------

* CSC, "Roihu" system description: https://docs.csc.fi/computing/systems-roihu/
* CSC, "Compiling on Roihu": https://docs.csc.fi/computing/compiling-roihu/
* CSC, "Available batch job partitions": https://docs.csc.fi/computing/running/batch-job-partitions/
* CSC, "Finland's National Supercomputer Roihu": https://csc.fi/en/our-expertise/high-performance-computing/roihu-supercomputer/
