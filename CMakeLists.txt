cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(recipe-09 LANGUAGES Fortran)

include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# add and build the boys function library
add_subdirectory(libboys)
# defines targets and sources
add_subdirectory(utility_library)

add_subdirectory(math_library)

add_subdirectory(src)

# contains an "external" library we will link to
# add_subdirectory(external)

# enable testing and define tests
enable_testing()
add_subdirectory(tests)

