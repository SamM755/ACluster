# Full-table end-to-end experiment

Rebuttal to Reviewer w6VF W3, Reviewer 1V54 W2/Q3, Reviewer 81fu W1 and Reviewer YZNw W3.

From the repository root:

```bash
python encoding_benchmark/scripts/run.py end-to-end
```

The experiment uses all rows and original columns of Musk, Crop, Gas, CT and BT for Chimp, Elf, Gorilla, RLE, Huffman, ALP, REGER, Vortex, KCluster and ACluster. Defaults are three warmups, ten measured repetitions and 10,000 rows per page.

ACluster and KCluster encode selected columns; ALP encodes the remaining columns in the same tuple order. Decoding preserves complete tuples and duplicate multiplicities. CR counts all serialized bytes, including metadata and binary64 corrections.

Encoding time includes column selection, preprocessing and complete stream generation, excluding CSV loading, file I/O and verification. Decoding time excludes verification. Output is written under `encoding_benchmark/output/`; `summary.csv` reports CR and mean encoding/decoding times in ns/row.

For Reviewer YZNw follow-up P3, use the per-run `SelectionMs` and `EncodeMs` fields. Group by dataset and method and compute `100 * mean(SelectionMs) / mean(EncodeMs)`.

Implementation: [pipeline.cpp](pipeline.cpp). Experiment driver: [experiments.cpp](experiments.cpp). Tests: [rebuttal index](../../RebuttalReproducibility.md#correctness-checks).
