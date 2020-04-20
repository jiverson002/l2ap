cmake_minimum_required(VERSION 3.0)
project(PackageTest C CXX)

include(${CMAKE_BINARY_DIR}/conan_paths.cmake)

find_package(L2AP CONFIG REQUIRED)

add_executable(example example.cpp)
target_link_libraries(example L2AP::l2ap OpenMP::OpenMP_CXX)
