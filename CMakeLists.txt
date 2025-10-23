cmake_minimum_required(VERSION 3.15)
set(PROJECT_NAME Helloworld_PETSc_MPI)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# for petsc
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# start a project
project(${PROJECT_NAME} LANGUAGES CXX C)

# Enable CTest for unit testing
enable_testing()

# Find required packages (PETSc should be discoverable via pkg-config)
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC IMPORTED_TARGET PETSc QUIET)
find_package(Kokkos REQUIRED)
find_package(MPI REQUIRED)

foreach(lib_name poisson_mms poisson_metrics poisson_solver)
    add_library(${lib_name} src/poisson/${lib_name}.cpp)
    target_include_directories(${lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
    if(NOT lib_name STREQUAL "poisson_mms")
        target_link_libraries(${lib_name} PUBLIC PkgConfig::PETSC MPI::MPI_CXX)
    endif()
endforeach()

# Poisson CLI test runner
add_executable(poisson_tests_cli src/poisson/tests/poisson_tests_cli.cpp)
target_link_libraries(poisson_tests_cli PRIVATE 
    poisson_mms poisson_metrics poisson_solver 
    PkgConfig::PETSC MPI::MPI_CXX)
target_include_directories(poisson_tests_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Add tests subdirectory
add_subdirectory(tests/unit/poisson)

