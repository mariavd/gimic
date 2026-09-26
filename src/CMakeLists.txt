include_directories(${CMAKE_BINARY_DIR})
include(ExternalProject)
add_subdirectory(libgimic)
add_subdirectory(fgimic)

configure_file(gimic.in gimic)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/gimic
    DESTINATION bin
    )

install(FILES fgimic/getkw.py
    DESTINATION ${PYTHON_SITE_INSTALL_DIR}
    )
