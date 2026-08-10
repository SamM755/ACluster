#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root/.deps/duckdb"
"$root/.build/test/unittest" '[cluster]'
"${CXX:-g++}" -std=c++17 -O1 -I"$root/include" -I"$root/codec" -Ithird_party/zstd/include \
  "$root/tests/codec_roundtrip.cpp" "$root/src/cluster_codec.cpp" "$root/src/cluster_compression.cpp" \
  "$root/codec/cluster_encoder_logic.cpp" "$root/src/cluster_hybrid.cpp" "$root/codec/alp.cpp" \
  "$root/.build/third_party/zstd/libduckdb_zstd.a" -o "$root/.build/codec_roundtrip"
"$root/.build/codec_roundtrip"
