import argparse
import csv
import ctypes
import os
import datetime
import json
from pathlib import Path
import shutil
import statistics

from benchmark_real import WORKSPACE, measure_queries, save, sha256, verify, check
from duckdb_local import Connection, DLL, ROOT, literal, version

FORMATS = ("raw_zstd", "parquet_zstd", "alp_zstd", "acluster_alp_zstd", "kcluster_alp_zstd")
LABELS = {"raw_zstd": "RAW+Zstd", "parquet_zstd": "Parquet+Zstd", "alp_zstd": "ALP+Zstd",
          "acluster_alp_zstd": "ACluster+ALP+Zstd", "kcluster_alp_zstd": "KCluster+ALP+Zstd"}

def summarize(report, directory):
    for item in report["datasets"]:
        for metrics in item["formats"].values():
            runs = metrics["runs"]
            if not runs:
                continue
            assert len({run["bytes"] for run in runs}) == 1
            assert len({run["artifact_sha256"] for run in runs}) == 1
            metrics.update(bytes=runs[0]["bytes"], compression_ratio=item["typed_bytes"] / runs[0]["bytes"],
                           write_ns=statistics.median(run["write_ns"] for run in runs),
                           scan_ns=statistics.median(ns for run in runs for ns in run["queries"]["full_scan"]["samples_ns"]),
                           bitwise_tuple_multiset="PASS")
    save(report, directory / "report.json")
    lines = ["# Full-column DuckDB experiment", "", f"Status: {report['status']}", "",
             "All original rows and columns are retained. Selected column IDs come from the existing manifest.",
             "Selected columns use ACluster/KCluster. Remaining columns use the repository's scalar ALP implementation in the same output row order.",
             "All file formats are written and queried through the locally built DuckDB C API.", "",
             "Each cell: complete file size (MB) / write time (s) / scan-and-aggregate time (ms).", "",
             "| Dataset (all/selected columns) | " + " | ".join(LABELS[f] for f in FORMATS) + " |",
             "|---|" + "---:|" * len(FORMATS)]
    for item in report["datasets"]:
        cells = []
        for name in FORMATS:
            m = item["formats"].get(name, {})
            cells.append(f"{m['bytes']/1e6:.3f} / {m['write_ns']/1e9:.3f} / {m['scan_ns']/1e6:.2f}" if "bytes" in m else "pending")
        lines.append(f"| {item['dataset']} ({item['dimension']}/{len(item['selected_columns'])}) | " + " | ".join(cells) + " |")
    if report["status"] == "PASS":
        cells = []
        for name in FORMATS:
            metrics = [item["formats"][name] for item in report["datasets"]]
            cells.append(f"{statistics.mean(m['bytes'] for m in metrics)/1e6:.3f} / "
                         f"{statistics.mean(m['write_ns'] for m in metrics)/1e9:.3f} / "
                         f"{statistics.mean(m['scan_ns'] for m in metrics)/1e6:.2f}")
        lines += ["| Avg. | " + " | ".join(cells) + " |", ""]
    lines += ["", "## Protocol", "", *[f"- {value}" for value in report["protocol"].values()], "",
              "The complete file accounts for column IDs, both codec streams, correction/exception data, checksums, framing and Zstd.",
              "ALP denotes the repository scalar implementation with local bit-exact and bounds-safety fixes, not DuckDB's native vectorized ALP storage codec.",
              "Avg. is an unweighted arithmetic mean of the displayed dataset metrics, not a throughput or geometric-mean ratio.", ""]
    (directory / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", default="Musk,Crop,Gas,CT,BT")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--cpu", type=int, default=None, help="Optional logical CPU affinity")
    args = parser.parse_args()
    assert args.repeats >= 1
    if args.cpu is not None:
        if args.cpu < 0:
            raise ValueError('CPU index must be nonnegative')
        if os.name == 'nt':
            kernel = ctypes.WinDLL('kernel32', use_last_error=True)
            kernel.GetCurrentProcess.restype = ctypes.c_void_p
            kernel.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
            if not kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(), 1 << args.cpu):
                raise ctypes.WinError(ctypes.get_last_error())
        elif hasattr(os, 'sched_setaffinity'):
            os.sched_setaffinity(0, {args.cpu})
        else:
            raise RuntimeError('CPU affinity is unavailable')
    manifest = WORKSPACE / "data/selected_columns.csv"
    entries = list(csv.DictReader(manifest.open(encoding="utf-8-sig")))
    for entry in entries:
        entry["SourcePath"] = "data/raw/" + entry["Dataset"] + ".csv"
        entry["OutputPath"] = "data/final/" + entry["Dataset"] + ".csv"
    names = args.datasets.split(",")
    assert set(names) <= {"Musk", "Crop", "Gas", "CT", "BT"}
    entries = [next(e for e in entries if e["Dataset"] == name) for name in names]
    check(names, 'raw')
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8)))
    out = ROOT / "output" / ("full_hybrid_" + now.strftime("%Y%m%d_%H%M%S"))
    out.mkdir(parents=True)
    report = {"status": "RUNNING", "scope": "all_original_columns", "output_directory": out.relative_to(ROOT).as_posix(),
              "started_at": now.isoformat(), "duckdb_version": version().decode(),
              "manifest": manifest.relative_to(WORKSPACE).as_posix(), "manifest_sha256": sha256(manifest),
              "parameters": {"repeats": args.repeats, "cpu_affinity": args.cpu, "page_size": 10000, "k": 100, "seed": 20260929,
                             "zstd_level": 3, "threads": 1},
              "protocol": {
                  "input": "All original data/raw rows and columns parsed as DOUBLE once per dataset. CSV import excluded from write timings.",
                  "selection": "Fixed selected-column IDs from selected_columns.csv. No column-selection search is timed.",
                  "write": f"Median of {args.repeats} SQL COPY timings from the same in-memory table, including preprocessing, both codecs, exceptions, serialization, Zstd and file I/O/synchronization. Format order rotates each round. No untimed write warm-up.",
                  "scan": f"All-column SUM(ABS(column)) with COUNT(*), one warm-up and three measured warm-cache queries per file. Median of {args.repeats * 3} samples. Filtered scans are also retained.",
                  "correctness": "Every saved file passes exact binary64 complete-tuple multiset comparison using bidirectional EXCEPT ALL, preserving duplicates and cross-column associations.",
                  "layout": "Pages of at most 10,000 rows. Remaining columns follow the Cluster output row order. The stream stores column IDs but no original-row permutation.",
                  "build": f"Same locally built DuckDB library and Release -O1 extension for all formats. Single thread, optional process CPU affinity={args.cpu}, no frequency lock. Parquet layout differs from Cluster pages.",
                  "ratio": "CR = 8 * actual row count * full column count / complete file bytes.",
                  "comparison": "Full-column COPY and scan measurements."},
              "implementation_sha256": {}, "datasets": []}
    snapshots = [Path(__file__), ROOT / "scripts/benchmark_real.py",
                 ROOT / "scripts/duckdb_local.py"]
    for folder, pattern in (("src", "*.cpp"), ("include", "*.hpp"), ("codec", "*.cpp"), ("codec", "*.h")):
        snapshots += sorted((ROOT / folder).glob(pattern))
    for path in [DLL, *snapshots]:
        report["implementation_sha256"][str(path.relative_to(ROOT))] = sha256(path)
    print("RESULTS " + str(out), flush=True)
    save(report, out / "report.json")
    for entry in entries:
        name = entry["Dataset"]
        source = WORKSPACE / entry["SourcePath"]
        dim = int(entry["OriginalDim"])
        selected = [int(c) for c in entry["SelectedColumns"].split(";")]
        names = [f"c{c}" for c in range(dim)]
        directory = out / name
        directory.mkdir()
        item = {"dataset": name, "input_path": source.relative_to(WORKSPACE).as_posix(), "input_sha256": sha256(source),
                "dimension": dim, "selected_columns": selected,
                "formats": {label: {"runs": []} for label in FORMATS}}
        report["datasets"].append(item)
        print("LOAD " + name, flush=True)
        with Connection() as db:
            db.execute("SET temp_directory=" + literal(directory / "temp"))
            types = "{" + ",".join(literal(c) + ":'DOUBLE'" for c in names) + "}"
            _, load_ns = db.execute("CREATE TABLE source AS SELECT * FROM read_csv(" + literal(source) +
                                    ",header=false,auto_detect=false,delim=',',parallel=false,columns=" + types + ")")
            count, _ = db.execute("SELECT count(*) n FROM source")
            n = int(count[0]["n"])
            assert n == int(entry["Rows"])
            item.update(rows=n, typed_bytes=8*n*dim, csv_load_ns=load_ns)
            bits = ",".join("cluster_bits(" + c + ") AS " + c for c in names)
            db.execute("CREATE TEMP VIEW expected_bits AS SELECT " + bits + " FROM source")
            median, _ = db.execute("SELECT quantile_disc(c0,0.5) threshold FROM source")
            predicate = "c0>" + median[0]["threshold"]
            aggregates = "count(*) n," + ",".join("sum(abs(" + c + ")) s" + str(i) for i, c in enumerate(names))
            expected = {label: db.execute("SELECT " + aggregates + " FROM source" + suffix)[0]
                        for label, suffix in (("full_scan", ""), ("filtered_scan", " WHERE " + predicate))}
            item.update(aggregates=aggregates, predicate=predicate, expected_aggregates=expected)
            for repeat in range(args.repeats):
                order = FORMATS[repeat % len(FORMATS):] + FORMATS[:repeat % len(FORMATS)]
                for label in order:
                    print(f"WRITE {name} {label} round={repeat+1} rows={n} dim={dim}", flush=True)
                    parquet = label == "parquet_zstd"
                    path = directory / (label + f"_r{repeat+1}" + (".parquet" if parquet else ".cluster"))
                    if parquet:
                        options = "FORMAT parquet,COMPRESSION zstd,COMPRESSION_LEVEL 3,ROW_GROUP_SIZE 10000"
                        view = "read_parquet(" + literal(path) + ")"
                    else:
                        method = label.removesuffix("_zstd")
                        options = "FORMAT cluster,METHOD " + literal(method) + ",K 100,PAGE_SIZE 10000,SEED 20260929,COMPRESSION zstd,COMPRESSION_LEVEL 3"
                        if method.endswith("_alp"):
                            options += ",CLUSTER_COLUMNS " + literal(",".join(map(str, selected)))
                        view = "read_cluster(" + literal(path) + ")"
                    _, write_ns = db.execute("COPY source TO " + literal(path) + " (" + options + ")")
                    metrics = {"round": repeat+1, "path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size,
                               "write_ns": write_ns, "artifact_sha256": sha256(path)}
                    if not parquet:
                        records, _ = db.execute("SELECT * FROM cluster_metadata(" + literal(path) + ")")
                        metrics["pages"] = [{k: v if k in ("mode", "fallback_reason", "compression") else int(v)
                                             for k, v in record.items()} for record in records]
                        assert sum(p["rows"] for p in metrics["pages"]) == n
                    verify(db, view, names)
                    metrics["bitwise_tuple_multiset"] = "PASS"
                    metrics["queries"] = measure_queries(db, view, aggregates, predicate, expected, n)
                    item["formats"][label]["runs"].append(metrics)
                    summarize(report, out)
                    print(f"PASS {name} {label} round={repeat+1} MB={metrics['bytes']/1e6:.3f} write_s={write_ns/1e9:.3f}", flush=True)
    report["status"] = "PASS"
    report["completed_at"] = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
    summarize(report, out)
    print("COMPLETE " + str(out / "REPORT.md"), flush=True)

if __name__ == "__main__":
    main()
