#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cmake -S "$root/.deps/duckdb" -B "$root/.build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE='-O1 -DNDEBUG' \
  -DBUILD_UNITTESTS=ON -DENABLE_JEMALLOC=OFF -DDISABLE_BUILTIN_HTTPLIB=ON \
  -DBUILD_EXTENSIONS=cluster
cmake --build "$root/.build" --target shell unittest duckdb --parallel "${JOBS:-4}"
