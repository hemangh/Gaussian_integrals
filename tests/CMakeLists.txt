add_executable(overlap_test test_Gaussian_overlap.f90)

target_link_libraries(overlap_test PRIVATE Gaussian_integrals_lib)

add_test(
  NAME
    test_overlap_integrals
  COMMAND
    $<TARGET_FILE:overlap_test>
  )

add_executable(simpleH2_test test_simpleH2.f90 reader_mod.f90)

target_link_libraries(simpleH2_test Gaussian_integrals_lib)
  
add_test(
  NAME
    test_simpleH2
  COMMAND
    $<TARGET_FILE:simpleH2_test>
    )

    add_executable(simpleN2_test test_primitivesN2.f90 reader_mod.f90)
  
    target_link_libraries(simpleN2_test Gaussian_integrals_lib)
    
  add_test(
      NAME
        test_simpleN2
      COMMAND
        $<TARGET_FILE:simpleN2_test>
  )

  add_executable(boys_test2 test_boys.f90)
    
  target_link_libraries(boys_test2 Gaussian_integrals_lib)

  add_test(
    NAME
      test_boys
    COMMAND
      $<TARGET_FILE:boys_test2>
)