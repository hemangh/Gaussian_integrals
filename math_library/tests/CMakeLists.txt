add_executable(boys_test2 test_boys.f90)
    
  target_link_libraries(boys_test2 math_lib)

  add_test(
    NAME
      test_boys
    COMMAND
      $<TARGET_FILE:boys_test2>
)

add_executable(cg_test test_clebsh_gordan_coef.f90)

target_link_libraries(cg_test PRIVATE math_lib)

add_test(
  NAME
    test_clebsh_gordan_coef
  COMMAND
    $<TARGET_FILE:cg_test>
  )

add_executable(sphHY_test test_sphHY.f90)

target_link_libraries(sphHY_test PRIVATE math_lib)

add_test(
  NAME
    test_SphHY_coef
  COMMAND
    $<TARGET_FILE:sphHY_test>
  )

add_executable(matrix_op_test test_matrix_operations.f90)

target_link_libraries(matrix_op_test PRIVATE math_lib)

add_test(
  NAME
    test_matrix_operations
  COMMAND
    $<TARGET_FILE:matrix_op_test>
  )

  add_executable(lebedev_test test_lebedev.f90)

  target_link_libraries(lebedev_test PRIVATE math_lib)
  
  add_test(
    NAME
    lebedev_test
    COMMAND
      $<TARGET_FILE:lebedev_test>
    )

    add_executable(rsph_integral test_integral_3rspH.f90)

    target_link_libraries(rsph_integral PRIVATE math_lib)
    
    add_test(
      NAME
      rsph_integral
      COMMAND
        $<TARGET_FILE:rsph_integral>
      )