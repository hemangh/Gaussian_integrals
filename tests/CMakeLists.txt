add_executable(fortran_test test_Gaussian_overlap.f90)

target_link_libraries(fortran_test Gaussian_integrals_lib)

add_test(
  NAME
    test_Gaussian_integrals_lib
  COMMAND
    $<TARGET_FILE:fortran_test>
  )