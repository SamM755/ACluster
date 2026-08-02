# ACluster

This repository contains the ACluster encoding benchmark, public datasets, and experiment drivers.

## Data

`encoding_benchmark/data/raw` contains the full-column public datasets. `encoding_benchmark/data/final` contains the columns used by the benchmark. The four native one-dimensional datasets are unchanged by column selection.

| Dataset | Raw dimensions | Final dimensions |
|---|---:|---:|
| SSD | 1 | 1 |
| BC | 1 | 1 |
| PI | 1 | 1 |
| WS | 1 | 1 |
| CT | 10 | 3 |
| Crop | 46 | 20 |
| Gas | 64 | 48 |
| BT | 77 | 49 |
| Musk | 166 | 111 |

Regenerate the final files with:

```bash
cd encoding_benchmark
python data/select_columns.py
```

The paper evaluates 12 datasets. This anonymous repository includes the nine publicly distributable datasets and omits TY, GW, and BX. The experiment entry points are unchanged, but aggregate results produced by this repository cover only the public subset.

The two BT files larger than 100 MiB are stored with Git LFS. After cloning, install Git LFS and retrieve them with `git lfs install` followed by `git lfs pull`.

## Build

Use MSYS2 MINGW64:

```bash
cd /d/path/to/ACluster/encoding_benchmark
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Run

Every experiment defaults to 100 measured runs and reports averages:

```bash
cd build
./benchmark.exe
```

Running without arguments is equivalent to `--experiment all`. A one-run smoke test can be selected without changing the source:

```bash
./benchmark.exe --experiment baseline --runs 1 --dataset SSD
```

Use `--data-root`, `--output-root`, and `--dataset` to override the defaults. Dataset names are case-sensitive. Results are written as CSV files under the output root.

## Experiment Entry Points

The paper mapping below follows the figure and table numbers in the compiled manuscript. Figures 6 and 7 are the primary encoding-baseline comparison.

| Entry | Command | Scope | Paper mapping | Main output |
|---|---|---|---|---|
| All experiments | `./benchmark.exe` | Runs every entry listed below | Orchestration only | `output/*` |
| Encoding baselines | `./benchmark.exe --experiment baseline` | Compression ratio and encoding/decoding time across Chimp, Elf, Gorilla, RLE, Huffman, ALP, REGER, Vortex, KCluster, and ACluster | Figures 6-7 (primary); Figures 9 and 13 use the same measurements | `output/baseline/results.csv` |
| Dimensionality | `./benchmark.exe --experiment dimensionality` | Evaluates increasing prefixes of the selected dimensions | Figure 8 | `output/dimensionality/dimensionality_raw_results.csv` |
| Order | `./benchmark.exe --experiment order` | Original order, random shuffling, and sorting by the first column | Table 2 | `output/order/results.csv` |
| One-dimensional special case | `./benchmark.exe --experiment rebuttal-1d` | Reproduces the SSD, BC, PI, and WS ACluster regression results | Figure 10 | `output/rebuttal_1d/results.csv` |
| Encoding components | `./benchmark.exe --experiment components` | Encoding-time and encoded-space breakdown for KCluster and ACluster | Figure 11 | `output/components/results.csv` |
| Outliers | `./benchmark.exe --experiment outlier` | ACluster under 0%, 0.1%, 0.5%, 1%, and 5% injected outliers | Table 4 | `output/outlier/results.csv` |
| KMeans comparison | `./benchmark.exe --experiment kmeans` | KCluster and ACluster compared with Euclidean KMeans encoding | Table 5 | `output/kmeans/results.csv` |
| Peak memory | `./benchmark.exe --experiment memory` | Peak encoding-memory comparison across all methods | Table 7 | `output/memory/memory_average.csv` |
| GZIP composition | `./benchmark.exe --experiment gzip` | Every supported encoding method followed by GZIP | Figure 14 | `output/gzip/results.csv` |
| Parameter sensitivity | `./benchmark.exe --experiment sensitivity` | KCluster `K`, page size, and pack size; ACluster page size and pack size | Not currently assigned a numbered paper result | `output/sensitivity/results.csv` |
| Reference comparison | `./benchmark.exe --experiment reference` | Adaptive ACluster references versus KCluster with fixed `K` values | Not currently assigned a numbered paper result | `output/reference/results.csv` |
| Subspace grouping | `./benchmark.exe --experiment subspace` | Full-space KCluster/ACluster versus grouped-dimension KClusterSub | Not currently assigned a numbered paper result | `output/subspace/subspace_group_raw_results.csv` |

The adjacent-cluster-scope experiment reported in Table 6 is intentionally not included because it requires changing ACluster's internal comparison scope.

## Common Examples

Run the primary paper experiment on all included public datasets:

```bash
./benchmark.exe --experiment baseline
```

Run the order experiment on one dataset:

```bash
./benchmark.exe --experiment order --dataset CT
```

Run the full suite once as a smoke test and write results to a separate directory:

```bash
./benchmark.exe --experiment all --runs 1 --output-root smoke_output
```
