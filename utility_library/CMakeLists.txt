add_library(utility_lib "")

target_sources(utility_lib
    PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/precision.f90
    ${CMAKE_CURRENT_LIST_DIR}/constants.f90
    ${CMAKE_CURRENT_LIST_DIR}/string_array_mod.f90
)