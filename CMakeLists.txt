cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME Helloworld_PETSc_MPI)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# for petsc
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)


# start a project
project(${PROJECT_NAME} LANGUAGES CXX C)

add_executable(ex1 src/main.cpp)


# find packages
include(FindPkgConfig REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(Kokkos REQUIRED)
find_package(MPI REQUIRED)


pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)

# link libraries
target_link_libraries(ex1 PkgConfig::PETSC Kokkos::kokkos MPI::MPI_CXX)
