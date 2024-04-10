add_library(math_lib "")

target_sources(math_lib
    PUBLIC
    conversions.f90
    special_functions.f90
    boys_function.f90
    factorials.f90
    lmindx_map.f90
    sphericalHarmonics.f90
    # clebsch_gordan.f
    clebsh_gordan.f90
    matrix_operations.f90
)

target_link_libraries(math_lib PUBLIC utility_lib -llapack -lblas)
enable_testing()
add_subdirectory(tests)