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

# 关掉msvc的一些warning
if (MSVC)
    add_compile_options(/wd4267) # size_t 转换的 warning
    add_compile_options(/wd4244) # double 转为 float 的 warning
    add_compile_options(/wd4305) # double 到 float 截断
    add_compile_options(/wd4819) # 编码相关的 warning
endif()

# find_package(OpenMP REQUIRED)
# if (OpenMP_FOUND)
#     message(STATUS "found openmp")
#     set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS}, ${OPENMP_C_FLAGS})
#     set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}, ${OPENMP_CXX_FLAGS})
# endif()

find_package(spdlog CONFIG REQUIRED)

message(STATUS "add submodule pmp-library")
add_subdirectory(external/pmp-library)
add_subdirectory(src)
