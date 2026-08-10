# Order heuristic

Rebuttal to Reviewer 81fu W2.

From the repository root:

```bash
python encoding_benchmark/scripts/run.py order
```

The runner uses all rows of BC and BT from `data/final`, with at most 10,000 rows per page and unchanged page membership. It compares ACluster with a cleanup pass that reassigns rows to lower-cost references and deletes a reference when its storage cost exceeds the added residual cost. Each reference is considered for removal at most once.

The default run uses one warmup and three measured repetitions per page. Encoding time includes preprocessing, reference selection, cleanup and serialization; CSV input, file writes and verification are excluded. Each page is decoded and checked as an exact binary64 tuple multiset.

Output `summary.csv` records reference counts, compressed bytes and encoding/decoding times. `SharedCR` uses a shared file header and page payloads for the rebuttal comparison; `IndependentCR` counts complete independent page streams.

To sort by the first column within each existing page:

```bash
python encoding_benchmark/scripts/run.py order --layout page-sort
```

The ACluster baseline rows from this command give the within-page sorting comparison.

Implementation: [cleanup.cpp](cleanup.cpp). Experiment driver: [experiments.cpp](experiments.cpp). Tests: [rebuttal index](../../RebuttalReproducibility.md#correctness-checks).
