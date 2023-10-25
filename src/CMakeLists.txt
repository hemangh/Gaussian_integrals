add_library(Gaussian_integrals_lib "")

target_sources(Gaussian_integrals_lib
    PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/precision.f90
    ${CMAKE_CURRENT_LIST_DIR}/constants.f90
    ${CMAKE_CURRENT_LIST_DIR}/contracted_Gaussian_type.f90
    ${CMAKE_CURRENT_LIST_DIR}/Gaussian_Overlap.f90
    ${CMAKE_CURRENT_LIST_DIR}/kinetic_integral.f90
    ${CMAKE_CURRENT_LIST_DIR}/boys_function.f90
    ${CMAKE_CURRENT_LIST_DIR}/nuclear_attraction_integral.f90
    ${CMAKE_CURRENT_LIST_DIR}/electron_repulsion_integral.f90
    )


