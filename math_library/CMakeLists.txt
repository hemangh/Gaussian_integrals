add_library(math_lib "")

target_sources(math_lib
    PUBLIC
    special_functions.f90
    boys_function.f90
    factorials.f90
    lmindx_map.f90
    sphericalHarmonics.f90
)

target_link_libraries(math_lib PUBLIC utility_lib)