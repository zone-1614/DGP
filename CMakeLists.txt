cmake_minimum_required (VERSION 3.5)

# set c++ standard
message(STATUS "set c++ standard")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set binary output paths
message(STATUS "set binary output path")
set(PROJECT_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

# set module path
message(STATUS "set module path")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

# set pmp-library options
set(PMP_BUILD_EXAMPLES OFF)
set(PMP_BUILD_TESTS OFF)
set(PMP_BUILD_DOCS OFF)

project ("DGP" LANGUAGES CXX)

# set windows specific definition :)
if(WIN32)
    add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX)
endif()

find_package(spdlog CONFIG REQUIRED)

message(STATUS "add submodule pmp-library")
add_subdirectory(external/pmp-library)
add_subdirectory(src)
