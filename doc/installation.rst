

Installation
============

Before compiling GIMIC you need to make sure that you have installed
ideally all of the packages collected in
the ``requirements.txt`` file.
You need minimum the ones listed below. ::

* cython
* numpy
* runtest == 2.3.2

A convenient way to install the packages listed in ``requirements.txt``
is to install the 
Anaconda2 (https://www.anaconda.com/distribution/) Python distribution
first. 
Then you can simply install all the packages listed in the file
``requirements.txt`` by using:: 

  $ conda install name-of-package

GIMIC requires CMake to configure and build. CMake is invoked via a front-end script called ``setup``::

  $ ./setup
  $ cd build
  $ make
  $ make install

To see all available options, run::

  $ ./setup --help

Branch "master"::
GIMIC requires CMake to configure and build.::

  $ mkdir build
  $ cd build
  $ cmake ../
  $ make
  $ make install

Test the installation with::

  $ cd build
  $ make test

Note, some tests may require Valgrind and will fail if this
debugger is not available. However, this is no need to worry if all
other tests pass. 


Parallelization
---------------

OpenMP parallelization is enabled by default. To build without it::

  $ ./setup --no-omp

At run time the number of threads is controlled by ``OMP_NUM_THREADS``.
BLAS is called from inside the OpenMP-parallel loops, so the ``gimic``
launcher sets ``OPENBLAS_NUM_THREADS``/``MKL_NUM_THREADS`` to 1 unless
they are already set; a threaded BLAS on top of OpenMP oversubscribes the
cores and runs several times slower.

MPI parallelization is in the works.


GPU offload
-----------

The evaluation of the current tensor (basis functions, density-matrix
products and contractions) can run on a GPU through OpenMP target
offload. With GNU compilers and the nvptx offload back end installed
(Debian/Ubuntu: ``gcc-<version>-offload-nvptx``)::

  $ ./setup --offload

For AMD GPUs pass ``--cmake-options="-DOFFLOAD_TARGET=amdgcn-amdhsa"``
as well. Other compilers need their offload flags in ``--extra-fc-flags``.

At run time the environment variable ``GIMIC_OFFLOAD`` selects the path:

* ``auto`` (default): use the GPU if one is present, otherwise the CPU code
* ``host``: run the offload code path on the CPU (for testing)
* ``off``: always use the CPU code

The GPU path currently requires a Cartesian basis (``spherical=off``, the
default) and prints the block size and device memory it uses. To check
a GPU build, ``tools/offload-check.sh`` runs the benzene 3D test case on
the CPU and on the GPU and compares the results.


Clusters with containerised conda (CSC Tykky)
---------------------------------------------

Some clusters provide conda only inside a container, for example CSC's
Tykky tool on Roihu (and formerly Puhti). ``conda-containerize`` builds
the environment into a container image and puts wrapper scripts for its
programs in ``<install dir>/bin``; running a wrapper starts the program
inside the container.

This affects GIMIC's ``gimic`` launcher, which is a Python script. When
GIMIC is built inside the container, CMake finds the container's own
Python, a path such as ``/CSC_TYKKY_xxxxxx/miniforge/envs/env1/bin/python3.12``,
and writes it into the first line of ``bin/gimic``. That path exists
only inside the container, so running ``gimic`` from the login shell or
a batch job fails. Pass the wrapper instead with ``--python-shebang``::

  $ ./setup --python-shebang=<install dir>/bin/python3

This sets the first line of ``bin/gimic`` and of the Python tools. If the
wrapper directory is always on ``PATH`` when GIMIC runs,
``--python-shebang="/usr/bin/env python3"`` works as well. (Older builds
without this option need the first line of ``build/bin/gimic`` edited by
hand to point to the wrapper; the edit is lost when GIMIC is rebuilt.)

Example files for the whole procedure are in ``container/tykky/``:
``gimic.yml`` (the conda environment, the same packages as
``requirements.txt``) and ``post-install.txt`` (clones and builds GIMIC
inside the container, in the installation directory). In
``post-install.txt``, replace ``INSTALL_DIR`` with the absolute path of the
installation directory, load the compiler, CMake and BLAS modules you want
to build with, then::

  $ module load tykky
  $ mkdir <install dir>
  $ conda-containerize new --prefix <install dir> container/tykky/gimic.yml
  $ conda-containerize update <install dir> --post-install container/tykky/post-install.txt
  $ export PATH="<install dir>/bin:$PATH"

GIMIC is then run as ``<install dir>/gimic/build/bin/gimic``. Check the
cluster's documentation for current module names; CSC recommends
installing into the project application directory (``/projappl``).


Installation on Stallo supercomputer
------------------------------------

::

  $ git clone git@github.com:qmcurrents/gimic.git
  $ cd gimic
  $ module load Python/2.7.12-foss-2016b
  $ module load CMake/3.7.2-foss-2016b
  $ virtualenv venv
  $ source venv/bin/activate
  $ pip install -r requirements.txt
  $ ./setup
  $ cd build
  $ make
  $ make install


Using BLAS routines
-------------------

The density-matrix contractions, which dominate the run time, use BLAS
(``dgemm``) when a BLAS library is found at configure time; this is the
default. Without one the code falls back to the ``matmul`` intrinsic,
which is several times slower, and ``./setup`` prints a warning. To
build without BLAS deliberately::

  $ ./setup --no-blas

With Intel compilers and MKL use::

  $ ./setup --fc=ifort --cc=icc --cxx=icpc --cmake-options="-D ENABLE_MKL_FLAG=ON"


Installation on a Mac
------------------------------

The main problem with installing GIMIC on a Mac is setting the
correct paths and force the program to find all libraries. 
Here we share some recommendations but 
can not guaratee that GIMIC will for sure work on your Mac. GIMIC has
been developed on a Ubuntu/Linux operating system. 
Note, there are for sure more elegant ways to get GIMIC installed
on a Mac. If you figure them out, please share. 

*   Install git and cmake. 
*   Install Anaconda2 then you should be able to install most of the
    recommended packages listed in "requirements.txt". 
*   For installing the C++ (gcc) and Fortran (gfortran) compilers you
    can use Xcode, Brew or MacPorts. You just need to make sure that
    GIMIC is able to find them. A way to solve this is editing the
    ".bashrc" file. 
*   Check if you have the following paths in your ".bashrc". If not
    add them. 
        * export PATH=/Users/your_username:/Users/your_username/anaconda/bin:$PATH 
        * PATH="/Applications/CMake.app/Contents/bin":"$PATH"
        * export PATH=/Users/your_username/gimic/build/bin:$PATH
        * export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"







