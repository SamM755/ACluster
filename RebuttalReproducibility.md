# Rebuttal reproduction

Run commands from the repository root after the [benchmark build](readme.md#build) or [DuckDB build](extension/duckdb/README.md#build-and-test).

## Experiments

| Rebuttal location | Script or implementation | Command | Output files |
|---|---|---|---|
| Reviewer w6VF W3; Reviewer 1V54 W2/Q3; Reviewer 81fu W1; Reviewer YZNw W3 | [End-to-end experiment](encoding_benchmark/end_to_end/experiments.cpp) | `python encoding_benchmark/scripts/run.py end-to-end` | `encoding_benchmark/output/end-to-end_*/summary.csv` |
| Reviewer 1V54 W1/Q1; Reviewer YZNw W2/Q3 and follow-up P2 | [Optimal experiment](encoding_benchmark/optimal/experiments.cpp) | `python encoding_benchmark/scripts/run.py optimal` | `encoding_benchmark/output/optimal_*/summary.csv` |
| Reviewer 81fu W2, cleanup heuristic | [Order experiment](encoding_benchmark/order/experiments.cpp), [cleanup](encoding_benchmark/order/cleanup.cpp) | `python encoding_benchmark/scripts/run.py order` | `encoding_benchmark/output/order_*/summary.csv` |
| Reviewer 81fu W2, within-page sorting | [Order experiment](encoding_benchmark/order/experiments.cpp) | `python encoding_benchmark/scripts/run.py order --layout page-sort` | `encoding_benchmark/output/order_*/summary.csv` |
| Reviewer YZNw W1/Q1, DuckDB table | [benchmark_full_hybrid.py](extension/duckdb/scripts/benchmark_full_hybrid.py) | `python extension/duckdb/scripts/benchmark_full_hybrid.py` | `extension/duckdb/output/full_hybrid_*/REPORT.md`, `report.json` |
| Reviewer YZNw follow-up P3, selection-cost table | [End-to-end experiment](encoding_benchmark/end_to_end/experiments.cpp) | `python encoding_benchmark/scripts/run.py end-to-end` | `encoding_benchmark/output/end-to-end_*/<dataset>/all_runs.csv`: `SelectionMs`, `EncodeMs` |

## Correctness checks

| Rebuttal location | Test code | Command | Output |
|---|---|---|---|
| Reviewer YZNw W1/Q1 and follow-up P3, binary64 | [codec_roundtrip.cpp](extension/duckdb/tests/codec_roundtrip.cpp), [cluster_roundtrip.test](extension/duckdb/tests/cluster_roundtrip.test) | `bash extension/duckdb/scripts/test.sh` | Console test summary |
| Reviewer YZNw follow-up P1/P3, complete tuples | [cluster_hybrid.test](extension/duckdb/tests/cluster_hybrid.test) | `bash extension/duckdb/scripts/test.sh` | Console test summary |
| Reviewer YZNw follow-up P1/P3, file round trips | [integration_test.py](extension/duckdb/scripts/integration_test.py) | `python extension/duckdb/scripts/integration_test.py` | `extension/duckdb/output/integration_*/report.json` |
| Reviewer 81fu W1; Reviewer YZNw W3 and follow-up P1 | [test_end_to_end_pipeline](encoding_benchmark/end_to_end/pipeline.cpp) | `ctest --test-dir encoding_benchmark/build -R table_roundtrip --output-on-failure` | Console test summary |
| Reviewer 1V54 W1/Q1; Reviewer YZNw W2/Q3 and follow-up P2 | [test_optimal_1d](encoding_benchmark/optimal/test.cpp) | `ctest --test-dir encoding_benchmark/build -R table_roundtrip --output-on-failure` | Console test summary |
| Reviewer 81fu W2 | [experimentSelfTest](encoding_benchmark/order/cleanup.cpp) | `ctest --test-dir encoding_benchmark/build -R order_cleanup --output-on-failure` | Console test summary |
