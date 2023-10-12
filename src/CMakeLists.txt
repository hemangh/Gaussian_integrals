add_library(Gaussian_integrals_lib "")

target_sources(Gaussian_integrals_lib
    PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/precision.f90
    ${CMAKE_CURRENT_LIST_DIR}/constants.f90
    ${CMAKE_CURRENT_LIST_DIR}/Gaussian_type.f90
    ${CMAKE_CURRENT_LIST_DIR}/Gaussian_Overlap.f90
    )


