cmake_minimum_required(VERSION 3.12)

project(CompressionBenchmark C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG") 
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_DEBUG "-g")

file(GLOB UTILS_SOURCES "utils/*.cpp" "utils/*.c")
file(GLOB ALGO_SOURCES "algorithms/*.cpp")

add_executable(
    benchmark
    main.cpp
    ${UTILS_SOURCES}
    ${ALGO_SOURCES}
)

target_include_directories(benchmark PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(benchmark PRIVATE stdc++fs) # GCC/Clang

message(STATUS "CMake configuration complete. You can now build the 'benchmark' target.")