include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})

set(sources
    globals.f90
    basis.f90
    bfeval.f90
    cao2sao.f90
    caos.f90
    dens.f90
    factorial.f90
    gaussint.f90
    gtodefs.f90
    intgrl.f90
    jtensor.F90
    jtensor_offload.F90
    kinds.f90
    lip.f90
    reorder.f90
    settings.f90
    teletype.f90
    tensor.f90
    timer.f90
    acid.f90
    )


add_library( gimic2
    ${sources}
    )
set_property( TARGET gimic2 PROPERTY POSITION_INDEPENDENT_CODE ON )

if(ENABLE_BLAS)
    target_compile_options(
        gimic2
        PRIVATE
            -DHAVE_BLAS
        )
endif()

option(ENABLE_MKL_FLAG "Enable MKL compile flag" OFF)
if(ENABLE_MKL_FLAG)
    target_compile_options(
        gimic2
        PRIVATE
            -DHAVE_BLAS
        )
endif()

install(TARGETS gimic2
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

