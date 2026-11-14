# The interactive scripts, configured for this installation and installed
# next to gimic: they are wrappers that call gimic-prepare by its full path.

set(GIMIC_PREPARE ${CMAKE_INSTALL_FULL_BINDIR}/gimic-prepare)

file(GLOB _script_templates RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.sh.in)
set(_scripts)
foreach(_t ${_script_templates})
    string(REGEX REPLACE "\\.in$" "" _s ${_t})
    configure_file(src/${_t} ${_s} @ONLY)
    list(APPEND _scripts ${CMAKE_CURRENT_BINARY_DIR}/${_s})
endforeach()

install(PROGRAMS ${_scripts}
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

# the aliases, for reference
install(FILES useful-aliases
    DESTINATION ${CMAKE_INSTALL_DATADIR}/gimic/jobscripts
    )
