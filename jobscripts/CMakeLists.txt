# The interactive scripts, configured for this installation and installed
# next to gimic: the profile wrappers call gimic-profile by its full path,
# 3D-run.sh reads its helper files from the data directory.

set(GIMIC_PROFILE ${CMAKE_INSTALL_FULL_BINDIR}/gimic-profile)
set(SCRIPTS_DIR ${CMAKE_INSTALL_FULL_DATADIR}/gimic/jobscripts)

file(GLOB _script_templates RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.sh.in)
set(_scripts)
# Only @SCRIPTS_DIR@ and @GIMIC_PROFILE@ are replaced: configure_file would
# also empty the placeholders that 3D-run.sh fills in at run time (@MINX@ ...).
foreach(_t ${_script_templates})
    string(REGEX REPLACE "\\.in$" "" _s ${_t})
    file(READ ${CMAKE_CURRENT_SOURCE_DIR}/src/${_t} _text)
    string(REPLACE "@SCRIPTS_DIR@" "${SCRIPTS_DIR}" _text "${_text}")
    string(REPLACE "@GIMIC_PROFILE@" "${GIMIC_PROFILE}" _text "${_text}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${_s} "${_text}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS src/${_t})
    list(APPEND _scripts ${CMAKE_CURRENT_BINARY_DIR}/${_s})
endforeach()

install(PROGRAMS ${_scripts}
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

# helper files of 3D-run.sh, and the aliases for reference
install(FILES src/functions-def src/gimic.cdens
    DESTINATION ${CMAKE_INSTALL_DATADIR}/gimic/jobscripts/src
    )
install(FILES useful-aliases
    DESTINATION ${CMAKE_INSTALL_DATADIR}/gimic/jobscripts
    )
