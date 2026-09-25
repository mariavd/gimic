# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE

#.rst:
#
# Find and link to BLAS.
#
# Variables defined::
#
#   BLAS_FOUND
#   BLAS_LIBRARIES
#   BLAS_INCLUDE_DIR
#
# autocmake.yml configuration::
#
#   docopt: "--blas Find and link to BLAS [default: False]."
#   define: "'-DENABLE_BLAS={0}'.format(arguments['--blas'])"

# GIMIC: BLAS is on by default and optional. Without it the density-matrix
# contractions fall back to the matmul intrinsic, which is much slower.

option(ENABLE_BLAS "Find and link to BLAS" ON)

if(ENABLE_BLAS)
    find_package(BLAS)
    if(NOT BLAS_FOUND)
        message(WARNING "No BLAS library found: building with the (much slower) matmul fallback. Install one (e.g. OpenBLAS) or configure with --no-blas to silence this warning.")
        set(ENABLE_BLAS OFF)
    endif()
endif()
