function(gimic_test _name)
  add_test(
    NAME
      ${_name}
    COMMAND
      ./test --binary-dir=${PROJECT_BINARY_DIR}/bin
             --work-dir=${PROJECT_BINARY_DIR}/test/${_name}
             --verbose
    WORKING_DIRECTORY
      ${PROJECT_SOURCE_DIR}/test/${_name}
    )
endfunction()


# Run a case with several OpenMP thread counts; pointwise output must be
# bitwise identical, integrals must agree to 1e-10. See test/threads/test.
function(gimic_thread_test _name _case _threads)
  add_test(
    NAME
      threads/${_name}
    COMMAND
      ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/threads/test
             --binary-dir=${PROJECT_BINARY_DIR}/bin
             --work-dir=${PROJECT_BINARY_DIR}/test/threads/${_name}
             --case=${PROJECT_SOURCE_DIR}/test/${_case}
             --threads=${_threads}
    )
endfunction()

gimic_test(benzene/integration-gauss)
gimic_test(benzene/integration-lobatto)
gimic_test(benzene/vectors)
gimic_test(benzene/2d)
gimic_test(benzene/3d)
gimic_test(benzene/keyword-magnet)
gimic_test(benzene/3d-keyword-magnet)
gimic_test(benzene/2d-keyword-magnet)
gimic_test(benzene/keyword-spacing)
gimic_test(benzene/keyword-rotation)
gimic_test(benzene/keyword-rotation_origin)
gimic_test(benzene/keyword-radius)
gimic_test(benzene/int-jmod-radius)
gimic_test(benzene/int-jmod-radius-excluded)
gimic_test(benzene/int-grid-bond-even)
gimic_test(benzene/int-grid-std-gauss)
gimic_test(benzene/profile)
gimic_test(benzene/profile-orientation)
gimic_test(benzene/profile-default-fixpoint)

# the gimic-prepare command end to end (needs the gimic launcher and python3)
add_test(
  NAME prepare-tool
  COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/prepare-tool/test
          --binary-dir=${PROJECT_BINARY_DIR}/bin
          --work-dir=${PROJECT_BINARY_DIR}/test/prepare-tool
  )
gimic_test(benzene/profile-squares)
gimic_test(benzene/int-cdens)
gimic_test(benzene/diamag-off)
gimic_test(benzene/paramag-off)
gimic_test(benzene/giao-test)
gimic_test(c4h4/integration)
gimic_test(c4h4/read-grid)
gimic_test(open-shell/3d)
gimic_test(open-shell/integration)
gimic_test(benzene/skip-jmod-integration)
gimic_test(benzene/magnetizability)
gimic_test(benzene/tropicity-cdens)
gimic_test(benzene/tropicity-capture)
gimic_test(benzene/tropicity-integration)
gimic_test(cot/tropicity-integration)
gimic_test(benzene/divj)
gimic_test(c4h4/divj-file)
gimic_test(open-shell/divj)
# calc=divj against itself and against calc=cdens; see test/divj-consistency/test
add_test(
  NAME
    divj-consistency
  COMMAND
    ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/divj-consistency/test
           --binary-dir=${PROJECT_BINARY_DIR}/bin
           --work-dir=${PROJECT_BINARY_DIR}/test/divj-consistency
           --case=${PROJECT_SOURCE_DIR}/test/benzene
  )
# binary .vti files against the ascii ones, through vti2csv; see test/vti-binary/test
add_test(
  NAME
    vti-binary
  COMMAND
    ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/vti-binary/test
           --binary-dir=${PROJECT_BINARY_DIR}/bin
           --work-dir=${PROJECT_BINARY_DIR}/test/vti-binary
           --case=${PROJECT_SOURCE_DIR}/test/benzene
  )
# a serial build would pass these trivially
if(ENABLE_OPENMP)
  gimic_thread_test(cdens-closed-shell benzene/3d 1,4,16)
  gimic_thread_test(cdens-open-shell open-shell/3d 1,4)
  gimic_thread_test(integration-closed-shell benzene/integration-gauss 1,4,16)
  gimic_thread_test(integration-open-shell open-shell/integration 1,4,16)
  gimic_thread_test(tropicity-cdens benzene/tropicity-cdens 1,4)
  gimic_thread_test(tropicity-capture benzene/tropicity-capture 1,4)
  gimic_thread_test(tropicity-integration benzene/tropicity-integration 1,4)
  gimic_thread_test(profile benzene/profile 1,4)
  gimic_thread_test(divj-closed-shell benzene/divj 1,4,16)
endif()
