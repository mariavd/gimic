# Current profiles: bugs found, port to Fortran, and a friendlier interface

Working notes and plan, October 2026. Everything marked *verified* was
reproduced in this session with a fresh build of GIMIC on the benzene test
system (`test/benzene`, 4 cores, gfortran 13, OpenBLAS).

## 1. Summary

* The current-profile scripts run one `gimic` process per slice (400 to 500
  processes for a typical profile). Numerically, a profile is nothing more than
  the per-block contributions of the composite Gauss quadrature that GIMIC
  already uses along the width of a bond grid. One run over the whole plane with
  `nslices` Gauss blocks reproduces the sum of the separate slices to the last
  printed digit (section 3.1), so the port needs no new integration code, only
  per-column bookkeeping in `integrate_current` and an output routine.
* The sign problem and the fix-point problem are two separate bugs and both are
  now understood and fixed in the scripts (section 2). The sign of the
  integrated current does not depend on the fix point at all. It depends on the
  order of the two bond atoms relative to the field, and the documentation told
  users to order them such that the result comes out negative. The fix point
  decides on which side of the bond the integration starts, and the script
  computed it from undefined variables, so it was anchored at the coordinate
  origin instead of at the bond.
* With 4 OpenMP threads and OpenBLAS left at its default thread count a GIMIC
  integration is 2.4 times *slower* than serial. Pinning OpenBLAS to one
  thread gives the expected 3.2 times speed-up. This has to be part of any
  parallelisation plan (section 3.5).

## 2. Bugs

### 2.1 The fix point was computed from undefined variables (verified)

`jobscripts/src/current-profile-header` computed the fixed point with an awk
cross product and then added `Ax` and `Ay`, which were never passed to awk
(`-v A1x=... -v A2x=...` are the names in scope). awk treats them as 0, so the
fixed point was `B x (A2 - A1)` placed at the *coordinate origin*, not at the
bond. Only `Fz` was anchored to the bond (`Az` is defined inside the awk
program). For a molecule centred on the origin this happens to put F somewhere
sensible, for anything else F lands on a random side of the bond or close to the
bond axis ("Basis vectors are linearly dependent"). The squares, vertical and
single-integration scripts used `calculateFixedCoord`, which rotates atom 2 by
-10 degrees *about the coordinate origin*: the same class of error.

What the fix point actually does, from `src/fgimic/grid.f90`
(`setup_bond_grid`, lines 165 to 260):

```
d   = A2 - A1                 normal of the integration plane (bond direction, basv3)
n   = (A1 - F) x (A2 - F)     normal of the plane (A1, A2, F); "up" axis of the grid is basv1 = -n
in  = n x d                   basv2: the side on which the first component of width is measured
origin = A1 + distance*d - out*basv2 - up*basv1
```

The plane is perpendicular to the bond whatever F is. F only chooses the side
of the bond on which the integration starts, and the rotation of the rectangle
about the bond. Consequently F must lie on the same side of the bond as the
point the profile starts from (the ring centre). The natural choice is the
centroid of the atoms the user already enters to define the start of the
integration, and that is what the scripts now use (`fixPointAndOrientation` in
`functions-def`). If the user typed a distance instead, F is the bond midpoint
displaced by 1 bohr along `B x d`.

### 2.2 The sign (verified)

`integrate_current` (`src/fgimic/integral.f90`, line 139) integrates
`J . basv3`, the flux along the bond direction A1 to A2. Nine runs on the
benzene reference plane (bond 1-2, width [-2.2, 5.0], the GIMIC convention is
that the diatropic ring current is positive, reference value +11.802 nA/T):

| case | magnet | bond | fix point | J (nA/T) | in-side (basv2) |
|---|---|---|---|---|---|
| reference | `magnet_axis=X` | [1,2] | atom 4 | **+11.802** | (-1,0,0) towards the ring |
| explicit B, as reference | `[0,0,-1]` | [1,2] | atom 4 | +11.802 | (-1,0,0) |
| B reversed | `[0,0,+1]` | [1,2] | atom 4 | **-11.802** | (-1,0,0) |
| atoms swapped | `[0,0,-1]` | [2,1] | atom 4 | **-11.802** | (-1,0,0) |
| F on the far side of the bond | `[0,0,-1]` | [1,2] | (5,0,0) | +10.815 (plane mirrored, sweeps the whole ring) | (+1,0,0) |
| `magnet_axis=X`, atoms swapped | `magnet_axis=X` | [2,1] | atom 4 | +11.802 | (-1,0,0) |
| `magnet_axis=X`, F far side | `magnet_axis=X` | [1,2] | (5,0,0) | -10.815 (mirrored) | (+1,0,0) |
| automatic F = M + B x d | `[0,0,-1]` | [1,2] | (1.277,0,0) | +11.802 | (-1,0,0) |
| automatic F, atoms swapped | `[0,0,-1]` | [2,1] | (1.277,0,0) | -11.802 | (-1,0,0) |

Rules that follow (and that the vector algebra confirms):

1. With an explicit `magnet` vector the sign of the result flips with the order
   of the bond atoms and with the sign of B. It does not depend on F. It is
   positive for a diatropic current exactly when `B . n > 0`, i.e. when B is
   parallel to the normal of the plane (A1, A2, F) with F on the ring side.
2. The side of F decides where the integration starts (the "in" side); a wrong
   side mirrors the plane, which is a different integral, not a sign change.
3. With `magnet_axis=X` the field is set equal to `n`, so the result is
   independent of the atom order and only the side of F matters.

Why "always flipped": `doc/scripts.rst` said "the atomic indices should be
given in counter-clockwise manner". With the default field (0, 0, -1) and the
usual view from +z, counter-clockwise order gives `B . n < 0`, hence a
negative diatropic current, every time. Driving the unmodified
`current-profile-local.sh` on benzene (verified):

| atom order | fix point produced | result |
|---|---|---|
| 2, 1 (the documented order) | (-2.63, 0, 0) | correct plane, **net -11.72 nA/T**, paratropic outside, diatropic inside |
| 1, 2 | (2.63, 0, 0), outside the bond line | plane mirrored, starts 2.6 bohr outside the molecule and crosses the whole ring; net -0.39 nA/T |

After the fix both orders give the same profile: net +11.81 nA/T, paratropic
inside the ring, diatropic outside, the ring centre at x = 0 and the bond at
x = 2.28 bohr (coarse 0.2 bohr slices; the reference value is 11.80).

The fix in the scripts: F = centroid of the start atoms; the script reverses the
order of the bond atoms in the GIMIC input (`bond=[a2,a1]` and
`distance = |d| - distance`) when `B . n < 0`, prints what it did, and records
it in `calculation.dat` (`gimic input: bond=[...] distance=... reversed=1`).
With `magnet_axis=X` nothing is reversed, since there the result is invariant.

### 2.3 Other script bugs fixed in this session

* `current-profile-local-submit` and `jobscript.IN`: the last
  `nslices mod ntasks` slices were started in the background and never waited
  for; `gradient.sh` then ran while they were still computing (missing or
  truncated last slices in `current_profile.dat`). Added the `wait`, and the
  local script now also skips slices that already finished, like the cluster
  job script.
* `current-profile-local-submit` ran `rm -rf *dat *eps GIMIC.*` in the
  directory the user started the script from, deleting any `.dat` file of
  theirs. The cleanup is now restricted to the profile directory.
  Same for the squares scripts and the SLURM job scripts.
* `$parallel=$($nsteps)` (five files) is not an assignment; the clamp of the
  number of tasks never worked. Now `parallel=$nsteps`.
* The `[X]` option (`magnet_axis=X`) used `sed '/s/magnet=/#magnet=/'`, which
  is a syntax error (`/s/` is an address), and in `current-profile-header` it
  ran before the input template had been copied. Fixed in the current and
  vertical profile scripts.
* Rotation origin: the answer was read into `enableRotationerigin` and tested
  as `enableRotationOrigin`, so the `rotation_origin` keyword could never be
  enabled. Fixed (`useRotationOrigin`).
* Vertical, squares and single-integration scripts: the bond centre was stored
  in `Bx By Bz` but `userInp` measures the start and end distances from
  `bx by bz`, so the "in" and "out" distances entered as atom indices were
  measured from the coordinate origin. Verified on benzene: start came out as
  0 instead of 2.28 bohr. Fixed by using the lowercase names.
* Vertical and squares scripts never substituted the field into the input
  template (`@MFx@` versus `@Bx@` in `gimic.Inp`), so GIMIC could not parse
  their input at all. Fixed.
* Squares script: the height check used `up - down` although both are entered
  as positive magnitudes, so the default up = down always triggered "height is
  0" and asked again. Fixed to `up + down`, like the width check.

Not fixed, for the record:

* `gimic.last.inp` is only removed when the plane preview is requested;
  otherwise it is counted as a slice by the submit scripts and `gradient.sh`.
* `nsteps` is truncated to an integer, so the last partial slice is dropped and
  the preview plane end lies one slice beyond the last real slice.
* `compareEqual` uses `=` instead of `==` in awk and assigns to a local, i.e.
  it is a no-op; `compareSmaller` likewise.
* `crit_pts.sh` parses the atom indices out of the directory name and breaks
  when the suffix contains digits.
* The hard-coded path `/homeappl/home/mariavd/scripts/normalise_plane` in the
  automatic field option.
* The three copies of the squares scripts at the top of `jobscripts/` duplicate
  `jobscripts/src/`.

### 2.4 Issues on the Fortran side relevant to the port

* `magnet.f90`, `check_field` (line 68): if `B . basv3 > 0` the field is
  silently reversed ("Left handed coordinate system, reversing magnetic
  field", printed as INFO). For planar molecules with B perpendicular to the
  bond this never triggers, but for a tilted bond it flips the sign of a
  profile depending on which way the bond happens to point. The orthogonality
  warning on line 80 (`abs(x)-1.d0 > 1.d-12`) can never fire for unit vectors.
  The port should decide the sign explicitly (section 3.3) and turn the silent
  flip into a warning.
* `integral.f90`: the MPI path is incomplete. `schedule` splits the columns
  over the ranks, but the partial sums are never reduced (the `collect_sum`
  calls are commented out with a TODO on lines 163 to 167), so with
  `mpirun -np N` rank 0 prints its own share as the total. Per-column arrays
  (section 3.3) make the reduction a single `MPI_Reduce`.
* `grid.f90`: `radius` was uninitialised for `std` grids, so integration over
  a std grid gave zero (fixed, see 3.6).
* `settings%magnet` is not normalised in `get_magnet`, so `magnet=[0,0,-2]`
  doubles the current. The scripts always write unit vectors; the port should
  normalise.
* OpenBLAS threads versus OpenMP threads, measured on a 648 x 180 grid
  (116 640 points, benzene):

  | OMP_NUM_THREADS | OPENBLAS_NUM_THREADS | wall time |
  |---|---|---|
  | 1 | 1 | 6.4 s |
  | 1 | 4 | 6.6 s |
  | 4 | 1 | **2.0 s** |
  | 4 | 4 (the default) | **15.4 s** |
  | 2 | 1 | 3.3 s |

  The `dgemm` calls in `ctensor_batch` are small (64 points per block), so
  BLAS threading buys nothing and oversubscription is severe. No thread
  pinning exists in this checkout or on any branch of `mariavd/gimic`
  (checked with `git grep` and `git log --all`); if it was done elsewhere
  it has not reached this fork. The launcher
  (`src/fgimic/gimic.py.in`) should set `OPENBLAS_NUM_THREADS=1` and
  `MKL_NUM_THREADS=1` unless the user set them, or the code should call
  `openblas_set_num_threads(1)` when built with OpenMP.

### 2.5 Should the default field be (0, 0, +1) instead of (0, 0, -1)?

No. Reversing B reverses J everywhere (`J(-B) = -J(B)`), so it flips the sign
of the flux through a plane with a fixed normal, but it does not change which
side of the bond the integration starts on and it does not change the value
of the diatropic current. The convention "diatropic positive" is a statement
about the sense of circulation *relative to B*, and the scripts now impose it
by choosing the bond direction from `B . n`, so the sign of B has become
irrelevant for the profile. There are two reasons to keep (0, 0, -1):

* It is what `magnet_axis=X` produces for the documented benzene example
  (bond 1-2, fixpoint 4): `n = (A1 - F) x (A2 - F) = (0, 0, -1)`. That is
  almost certainly where the habit comes from: people copied the field GIMIC
  printed for that example into their explicit inputs. With that F the
  explicit (0, 0, -1) and `magnet_axis=X` are the same calculation.
* All `cdens` visualisations and the viewing rule in `doc/interpretation.rst`
  ("look along the field, diatropic is clockwise") assume it. Changing the
  default would mirror every new vector plot with respect to the published
  ones.

With the old scripts, (0, 0, +1) together with the documented counter-clockwise
atom order would indeed have given positive results, which is probably why
workflows that "always worked" differ from ones that "always came out
flipped". Both were relying on an accident. What does matter for non-planar
molecules is that B is perpendicular to the ring plane; its sign does not.

## 3. Porting the profile to Fortran

### 3.1 Why it is cheap: slices are Gauss blocks

`setup_gauss_data` (`src/libgimic/gaussint.f90`, line 267) splits each grid
axis of length `l` into `npts/gauss_order` blocks of equal length and places a
`gauss_order`-point Gauss rule in each block. The scripts choose
`spacing` such that every slice contains exactly one block of 9 points across
its width (0.02 bohr) and 648 points along the height. A single grid over the
whole plane with `npts(2) = 9 * nslices` has block boundaries at exactly the
slice boundaries and the same nodes and weights. Verified on benzene with 20
slices of 0.1 bohr (height 16 bohr, 648 points):

| | net | dia | para |
|---|---|---|---|
| sum of 20 separate `gimic` runs | -3.15353 | 1.70900 | -4.86252 |
| one run, grid 648 x 180 | -3.15355 | 1.70898 | -4.86252 |

The difference is the 5-decimal rounding of the printed slice values. Cost: the
20 separate runs took 16.5 CPU-seconds, the single run 7.0 s on one thread.
The gap is the per-process start-up (reading `XDENS` and the basis, Python
launcher), which grows with the molecule while the integration grows with the
grid, so the gain is larger for real systems. A production-like slab of
648 x 900 points (100 slices of 0.02 bohr) takes 32 s on one thread for
benzene, so a full 500-slice profile is about 3 minutes serial and under one
minute on 4 threads, in one process, with one `grid.xyz` of the whole plane.

The profile values are therefore available inside `integrate_current` for
free: the quantity the scripts extract from each slice output is
`sum over the slice's columns j of w_j * (sum over i of w_i J.n)`, and the
positive/negative split is per grid point, so it is additive over blocks too.

### 3.2 Input

Keep `calc=integral` and the `Grid(bond)` section, add two keywords and make
`fixpoint`/`fixcoord` optional:

```
Grid(bond) {
    type=gauss
    bond=[1,2]
    towards=[1,2,3,4,5,6]   # NEW: atoms whose centroid is the reference centre C (the ring).
                            #      Replaces fixpoint/fixcoord: F = C, and the sign convention
                            #      "diatropic positive" is applied automatically (3.3).
    distance=1.3147
    width=[-2.63, 8.0]      # unchanged: in (negative, towards C) and out
    height=[-8.0, 8.0]
    gauss_order=9
    slices=[0.02, 0.0]      # NEW: slice thickness along width and height (bohr); 0 = no slicing.
                            #      Sets npts = gauss_order * nint(length/slice) on that axis.
    spacing=[0.025, 0.0, 0.0]   # as today for the axes that are not sliced (here the height)
}
```

* `slices=[dw, 0]` gives the current profile along the plane,
  `slices=[0, dh]` the vertical profile, `slices=[dw, dh]` the squares map.
  The same accumulation serves all three (`squares-profile` and
  `vertical-profile` become one keyword each instead of 350-line scripts).
* `towards` can also take a point: `towards_coord=[x,y,z]`. `fixpoint` and
  `fixcoord` stay for backward compatibility and for deliberately tilted
  planes; with `towards` set, F = C.
* The Python side (`src/gimic.in`, `parse_keywords`/`check_grid`, lines 82
  to 102 and 188) needs the two keywords and the relaxed "fixpoint or
  fixcoord" check.
* `spacing` on a sliced axis is ignored with a note; `grid_points` on a sliced
  axis is an error.

### 3.3 Code changes

`src/fgimic/grid.f90`, `setup_bond_grid`:

* read `towards`; compute C; set `this%origin` (the fix point) to C before the
  existing geometry code. Keep everything else. Store C and a flag in `grid_t`.
* when `slices` is set, compute `npts` here (instead of `setup_gauss_grid`
  rounding `l/spacing` up to a multiple of the order) so that the block edges
  are exactly where the user expects, and store the block count per axis.

`src/fgimic/integral.f90`, `integrate_current`:

* sign convention: with `towards`, define the reporting normal as
  `sign(B . n) * basv3` where `n` is the plane normal before the magnet is
  substituted into `ortho` (store it separately in `grid_t`; today `ortho` is
  overwritten by `magnet` on line 253). This is exactly the rule in 2.2 and is
  what the scripts now do by swapping the atoms; doing it in Fortran removes
  the swap, the `distance` mirroring, and the atom-order rule from the user's
  view. Print one line: "Current reported positive for flow in the diatropic
  sense around the reference centre". With `magnet_axis=X` the sign is already
  invariant and no change is needed.
* per-column accumulation: allocate `col(3, p2)` (net, positive, negative
  without the `w_j` weight) outside the parallel region and write
  `col(:, j)` inside the `!$OMP DO` loop over `j` (lines 120 to 152). Each
  column is owned by exactly one thread, so no reduction is needed; the
  totals become `sum_j w_j col(:, j)` after the loop and the existing
  `REDUCTION` clause can go. For the squares map, accumulate per
  (height block, column): `blk(3, nblk1, p2)`, also race-free. Memory is a
  few MB.
* progress: print every N columns with an estimate, as `jfield_eta` does for
  `cdens` (commit b29362d), because one profile run is now minutes to hours
  instead of hundreds of short jobs.
* `schedule(dynamic, 16)` on the column loop: columns far from the molecule
  are cheap because of screening, static scheduling leaves threads idle.

`src/fgimic/profile.f90` (new, small): given `grid` and `col`/`blk`, write

* `current_profile.dat`: the four columns the plotting and `crit_pts.sh`
  tools read today (distance from the "in" edge in bohr, net, dia, para in
  nA/T per slice), ordered from the in side to the out side. Note that in the
  grid frame `t = 0` is the *out* edge (origin = centre - out*basv2), so the
  slice index is reversed when writing.
* `profile.dat`: a header with the geometry (A1, A2, C, B, plane corners, the
  reporting normal), then per slice: start, centre and end position, net, dia,
  para in au and nA/T, and dJ/dx at the slice centre; plus the per-node values
  `t_j, col(:, j)` so that any sub-range can be re-integrated without rerunning.
* `squares_profile.dat`: the matrix `gradient-square.sh` plots, plus an
  `x y net dia para` long format.
* with `dryrun=on`: write the slice table and `grid.xyz` of the whole plane
  (the four corners replace `grid.check.xyz`).

`src/fgimic/parallel.F90`: replace the commented `collect_sum` calls by one
`MPI_Reduce` of `col` (and `blk`) onto rank 0. Each rank fills its `[lo, hi]`
columns and leaves the rest zero. This also fixes the broken MPI total of the
plain integral.

`src/fgimic/gimic.py.in`: pin BLAS threads (2.4).

### 3.4 Validation

* ctest `benzene/profile`: `slices=[0.1,0]` over width [-2.0, 0.0], 648 x 180,
  reference `current_profile.dat` from the 20 separate runs in 3.1 (they agree
  to 2e-5 nA/T) and the total from the existing `integration-gauss` test.
* ctest `benzene/profile-orientation`: the four combinations of atom order and
  field sign from 2.2 must give identical output with `towards=[1..6]`.
* ctest `benzene/squares`: `slices=[0.5,0.5]`, compare the matrix sum with the
  total.
* An MPI run (`GIMIC_LAUNCHER="mpirun -np 2"`) must reproduce the serial
  profile; today it does not even for the total.

### 3.5 Parallel execution after the port

* Default: one process, OpenMP over columns, BLAS pinned to one thread.
  Measured 3.2x on 4 cores; the column loop is embarrassingly parallel so this
  should hold to a full node.
* Clusters: one `sbatch` job with `--cpus-per-task=N` and
  `OMP_NUM_THREADS=N` replaces the job script that juggles hundreds of
  processes; several profiles become a job array. Multi-node only via MPI
  (after 3.3), which is rarely needed: a 500-slice benzene profile is a minute
  on 4 cores.
* Restart: no longer necessary for one run of minutes; for very large systems
  write `profile.dat` incrementally (columns are finished in order with
  `schedule(dynamic)` only approximately, so write completed blocks when all
  their columns are done) and add `--resume` later if anyone asks.

### 3.6 Planes in arbitrary positions, given by their corners

The request: define the integration rectangle directly by points in space,
e.g. a corner at the origin, "up" to (0, 0, h) and "out" to (l, 0, 0), with no
bond, fix point or rotation involved. GIMIC's `Grid(std)` is exactly this
description: `origin` is the corner P0, `ivec` and `jvec` are the directions
of the two edges, `lengths` their lengths, and `type=gauss` is accepted for
it. `integrate_current` works on any grid, with the normal `ivec x jvec`.
It did not work, though: `setup_bond_grid` initialises the cut-off `radius`
to 1e10 but `setup_std_grid` never set it, so for a std grid every point was
outside the (undefined) radius and the integral was silently **zero**. Fixed
in `new_grid` (`grid.f90`, default `radius = 1.d10` for all grid types) and
covered by the new test `benzene/int-grid-std-gauss`: the reference plane
rebuilt from its corners

```
Grid(std) {
    type=gauss
    origin=[7.277078, 0.0, -5.0]    # corner P0
    ivec=[0.0, 0.0, 1.0]            # P0 -> P1: the "up" edge, along B
    jvec=[-1.0, 0.0, 0.0]           # P0 -> P2: the "out" edge, along the profile
    lengths=[10.0, 7.2, 0.0]        # |P1 - P0|, |P2 - P0|
    gauss_order=9
    grid_points=[36, 36, 0]
}
```

gives 11.802083 nA/T (bond grid: 11.802134; the 5e-5 is the rounded corner
coordinates). Listing the corners in the other order (P2 on the other side)
gives -11.802083: for a std grid the sign is the flux along `ivec x jvec`,
nothing else, which is also the clean statement of the sign problem.

For the port this means the profile machinery should be built on the grid
axes, not on the bond grid:

* `slices` applies to axis 2 (and 1) of any grid type.
* `towards=[atoms]` or `towards_coord=[x,y,z]` applies to any grid type and
  defines the reporting normal as `sign((B x r) . basv3) * basv3`, with `r`
  from the reference centre C to the plane centre. For a bond grid with F on
  the ring side this is identical to the `B . n` rule of 2.2 (because
  `n` is parallel to `r x d` there); for a std grid it is the only way to get
  "diatropic positive" independent of the corner order.
* Convenience syntax in the parser: `corners=[P0, P1, P2]` (nine numbers) for
  `Grid(std)`, from which origin, ivec, jvec and lengths follow; and
  `from`/`to` points with `height` for a plane that contains a chosen line
  and B. Both are sugar on top of `Grid(std)`; `Grid(bond)` stays for the
  bond-centred workflow.
* The existing `coord1`/`coord2` route (`coord1 = P`, `coord2 = P + n`,
  `distance = 0`, verified to reproduce the reference plane) remains the way
  to ask for "a plane through P with normal n" within `Grid(bond)`.

### 3.7 Phases

1. `gimic.py.in` BLAS pinning; `check_field` flip becomes a warning. Half a
   day, independent of the rest.
2. `towards` + sign convention in `grid.f90`/`integral.f90`/`gimic.in`, with
   the orientation test. One to two days. This alone retires the fix-point
   question in every script.
3. `slices` + per-column accumulation + `profile.f90` + tests. Two to three
   days. The bash scripts then write one input file and call `gimic` once;
   `gradient.sh` becomes unnecessary (GIMIC writes `current_profile.dat`),
   `crit_pts.sh` and the plotting keep working on the same file.
4. MPI reduction and incremental output. One day, optional.

## 4. A friendlier interface

The interactive scripts ask about fifteen questions in a fixed order, most of
them "Do you accept X? Press [n] to change", and any mistake means starting
over (there is a partial retry menu for seven of the parameters). Concrete
suggestions, roughly in order of value for effort:

1. **One command instead of five scripts, in Python.** GIMIC already requires
   Python for the `gimic` launcher and the converters. A `gimic-profile`
   command with `argparse` gives non-interactive use for free:

   ```
   gimic-profile --bond 1 2 --ring 1 2 3 4 5 6 --out 8 --height 8 --slice 0.02 --field 0 0 -1
   gimic-profile --squares 0.25 ...
   gimic-profile --from current_profile_1.2/profile.json   # rerun or tweak
   ```

   Interactive prompting only for arguments that are missing, with the default
   shown in brackets and Enter accepting it (`in [2.78]:`), instead of
   "press n to change". Validation re-prompts instead of exiting
   ("Variable empty, please restart the script").
2. **Read XYZ directly.** The scripts need a Turbomole `coord` file and
   `coord.xyz`; GIMIC itself needs neither. Read the geometry from `MOL` (it is
   there already, in bohr) or from an XYZ file, accept Angstrom or bohr, and
   echo atoms as `C5 (2.277, 1.315, 0.000)` when the user picks them, so that
   a mistyped index is caught immediately.
3. **Remove the fix-point question.** With `towards` (or the script's new
   automatic F) the user never needs to know what a fix point is; keep an
   `--fixpoint` expert option for tilted planes.
4. **Automatic field direction without external binaries.** Fit a plane
   through the ring atoms (or all atoms) with a 3 x 3 SVD in numpy and use its
   normal, replacing the `maximise_projection`/`normalise_plane` programs and
   the hard-coded path. Offer `--field auto|z|x y z`.
5. **Preview before computing.** Write `plane.xyz` with the molecule, the four
   corners of the whole plane and the field arrow (one `gimic --dryrun` of the
   whole-plane input gives this), print the grid size and a run-time estimate,
   then ask once: run now, write a SLURM script, or quit. The current script
   asks the user to background the script with Ctrl-Z to look at the plane.
6. **One record of what was done.** Replace `calculation.dat` (free text,
   parsed with `sed` by four other scripts) with `profile.json`: inputs, the
   derived fix point and orientation, the GIMIC input file name, the GIMIC
   version. The plotting and critical-point tools read it, and
   `--from profile.json` reproduces a run.
7. **Plotting without gnuplot/ImageMagick/epstopdf.** `matplotlib` is
   available wherever numpy is; produce PNG and PDF with the ring centre at
   x = 0, vertical lines at the bond and at user-chosen atom centroids (what
   `plot-current-profile.sh` does interactively), and the integrated values
   between the zero crossings printed in the legend (what `anprofile` and
   `crit_pts.sh` do). Keep `current_profile.dat` so existing gnuplot
   workflows still work.
8. **Cluster use.** `--slurm` writes a job file next to the input with the
   partition, account and modules taken from a small user config
   (`~/.gimic-profile.toml`) and prints the `sbatch` line (`--submit` runs
   it). No editing of `jobscript-header` before `setup.sh`, and no `setup.sh`
   at all: the command is installed by CMake like `gimic` is.
9. **Interactive scripts as thin wrappers.** `gcurrent`, `gsquares` and the
   other aliases can call the Python command so that nobody has to relearn
   names.

Until the port is done, the bash scripts as fixed in this session are correct
for the common case (planar ring, start defined by atoms).

## 5. State of the repository after this session

Changed (uncommitted): `jobscripts/src/functions-def`,
`current-profile-header`, `current-profile-local-submit`,
`current-profile-cluster-submit`, `current-profile-cluster-restart.sh.in`,
`jobscript.IN`, `jobscript-squares.IN`, `vertical-profile-header`,
`squares-profile-header`, `squares-profile-local-submit`,
`squares-profile-cluster-submit`, `gimic-run.sh.in`, `gimic.Inp`,
`doc/scripts.rst`, and this file.

Verified: the current profile end to end on benzene for both atom orders
(local submit path, 52 slices, 4 tasks); the vertical, squares and single
integration scripts through their dry run and generated input; the single
integration run to completion (+11.42 nA/T, width [-2.28, 8], height +-10,
27 height points). Not verified: the SLURM paths (no scheduler here), the
plotting (no gnuplot here), the squares and vertical calculations beyond the
dry run, and the MPI statements, which are from reading the code.
