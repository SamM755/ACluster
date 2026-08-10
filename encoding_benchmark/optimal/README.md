# Optimal comparison

Rebuttal to Reviewer 1V54 W1/Q1 and Reviewer YZNw W2/Q3 and follow-up P2.

From the repository root:

```bash
python encoding_benchmark/scripts/run.py optimal
```

The experiment compares KOptimal, KCluster, AOptimal and ACluster on all rows of SSD, BC and PI. It uses independent pages of at most 10,000 rows, K=100, pack size 10, block size 10 and ten measured repetitions without warmups. K is capped at the distinct-value count.

CR includes complete page streams, headers and padding. Encoding time spans CSV-page reading through stream writing and closing, excluding verification and logging. Decoding time excludes verification. Each page is decoded and checked as an exact binary64 multiset.

[optimal_1d.cpp](optimal_1d.cpp) implements KOptimal and AOptimal; [cluster_1d.cpp](cluster_1d.cpp) implements KCluster and ACluster for this comparison. [experiments.cpp](experiments.cpp) runs the experiment. Output includes `runs.csv`, `pages.csv`, page streams and `summary.csv` under `encoding_benchmark/output/`. Times are in ns/row.

Correctness tests: [test_optimal_1d](test.cpp), run through the [rebuttal test command](../../RebuttalReproducibility.md#correctness-checks).
