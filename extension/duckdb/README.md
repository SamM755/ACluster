# DuckDB extension

Rebuttal to Reviewer YZNw W1/Q1 and follow-up P1/P3. The extension provides the `cluster` COPY format and `read_cluster` table function. The DuckDB commit is pinned in `dependencies.json`.

## Build and test

Requirements: Python 3.10+, Git, CMake, Ninja and a C++17 compiler. Run from the repository root in MSYS2 MINGW64; the build uses GCC Release `-O1`:

```bash
python extension/duckdb/scripts/prepare.py
bash extension/duckdb/scripts/build.sh
bash extension/duckdb/scripts/test.sh
python extension/duckdb/scripts/integration_test.py
```

`prepare.py` obtains the pinned DuckDB source and copies the extension and tests into it. Source and build files are kept in `.deps/` and `.build/`. The Python runners use the built DuckDB shared library; keep the compiler runtime libraries on `PATH`.

The [rebuttal test index](../../RebuttalReproducibility.md#correctness-checks) points to the binary64, SQL and complete-tuple tests.

## Reproduce the experiment

```bash
git lfs pull
python encoding_benchmark/scripts/check_data.py
python extension/duckdb/scripts/benchmark_full_hybrid.py
```

The experiment retains every row and original column of Musk, Crop, Gas, CT and BT. It compares RAW+Zstd, Parquet+Zstd, scalar ALP+Zstd, ACluster+ALP+Zstd and KCluster+ALP+Zstd. ACluster/KCluster encode selected columns, and ALP encodes the remaining columns in the same tuple order.

Complete file size includes metadata, corrections and framing. Write time includes preprocessing and I/O, excluding CSV import and selection search. Scans compute COUNT and SUM(ABS(column)) over all columns with one thread. The runner uses three writes per format and three measured scans per file after a warmup, reporting median write and scan times. The Avg row is the arithmetic mean across datasets.

Each file is checked as a binary64 multiset of complete tuples, preserving duplicate multiplicities. Precision mismatches use original-bit corrections; unsafe scaling and non-finite values use raw storage.

Results are written under `extension/duckdb/output/`. Column-selection timing for follow-up P3 is measured by the [full-table end-to-end experiment](../../encoding_benchmark/end_to_end/README.md).
