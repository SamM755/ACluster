cmake_minimum_required(VERSION 3.12)

project(CompressionBenchmark C CXX)

enable_testing()
file(GLOB TABLE_CODEC_SOURCES CONFIGURE_DEPENDS "end_to_end/codec/algorithms/*.cpp" "end_to_end/codec/utils/*.cpp")
add_executable(table_benchmark end_to_end/main.cpp end_to_end/io.cpp end_to_end/pipeline.cpp
    end_to_end/experiments.cpp optimal/optimal_1d.cpp optimal/cluster_1d.cpp
    optimal/pipeline.cpp optimal/experiments.cpp optimal/test.cpp ${TABLE_CODEC_SOURCES})
target_include_directories(table_benchmark PRIVATE end_to_end end_to_end/codec end_to_end/codec/algorithms optimal)
target_compile_features(table_benchmark PRIVATE cxx_std_17)
add_test(NAME table_roundtrip COMMAND table_benchmark --experiment test)
add_executable(order_benchmark order/main.cpp order/experiments.cpp order/cleanup.cpp algorithms/cluster_encoder_logic.cpp)
target_include_directories(order_benchmark PRIVATE . order)
target_compile_definitions(order_benchmark PRIVATE ACLUSTER_BUFFER_EXPERIMENT)
target_compile_features(order_benchmark PRIVATE cxx_std_17)

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

add_test(NAME order_cleanup COMMAND order_benchmark self-test)
