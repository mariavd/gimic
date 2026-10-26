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
# a serial build would pass these trivially
if(ENABLE_OPENMP)
  gimic_thread_test(cdens-closed-shell benzene/3d 1,4,16)
  gimic_thread_test(cdens-open-shell open-shell/3d 1,4)
  gimic_thread_test(integration-closed-shell benzene/integration-gauss 1,4,16)
  gimic_thread_test(integration-open-shell open-shell/integration 1,4,16)
endif()
