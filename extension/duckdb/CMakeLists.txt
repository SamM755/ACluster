cmake_minimum_required(VERSION 3.14)
project(ClusterExtension)
include_directories(include codec)
build_extension_library(cluster NO_LOADABLE
  src/cluster_extension.cpp src/cluster_codec.cpp src/cluster_hybrid.cpp src/cluster_compression.cpp
  codec/cluster_encoder_logic.cpp codec/alp.cpp)
target_compile_features(cluster_extension PRIVATE cxx_std_17)
target_link_libraries(cluster_extension duckdb_zstd)
