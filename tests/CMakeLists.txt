add_executable(fortran_test test_Gaussian_overlap.f90)

target_link_libraries(fortran_test Gaussian_integrals_lib)

add_executable(simpleH2_test test_simpleH2.f90)

target_link_libraries(simpleH2_test Gaussian_integrals_lib)

add_test(
  NAME
    test_Gaussian_integrals_lib
  COMMAND
    $<TARGET_FILE:fortran_test>
  )
  add_test(
  NAME
    test_simpleH2
  COMMAND
    $<TARGET_FILE:simpleH2_test>
    )
