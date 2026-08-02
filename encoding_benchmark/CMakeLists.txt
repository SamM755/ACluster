cmake_minimum_required(VERSION 3.12)

project(CompressionBenchmark C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_DEBUG "-g")

file(GLOB UTILS_SOURCES CONFIGURE_DEPENDS "utils/*.cpp" "utils/*.c")
file(GLOB ALGO_SOURCES CONFIGURE_DEPENDS "algorithms/*.cpp")
file(GLOB EXPERIMENT_SOURCES CONFIGURE_DEPENDS "experiments/*.cpp")

add_executable(
    benchmark
    main.cpp
    benchmark_core.cpp
    ${EXPERIMENT_SOURCES}
    ${UTILS_SOURCES}
    ${ALGO_SOURCES}
)

target_include_directories(benchmark PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(benchmark PRIVATE stdc++fs)
if(WIN32)
    target_link_libraries(benchmark PRIVATE psapi)
endif()
