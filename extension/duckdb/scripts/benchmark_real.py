import argparse
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

from duckdb_local import Connection, DLL, ROOT, literal, version

WORKSPACE = ROOT.parents[1] / "encoding_benchmark"
sys.path.insert(0, str(WORKSPACE / 'scripts'))
from check_data import check

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def check_aggregates(expected, actual):
    assert actual and actual[0].keys() == expected[0].keys(), "Missing aggregate result"
    for key in expected[0]:
        lhs, rhs = expected[0][key], actual[0][key]
        if lhs is None or rhs is None:
            assert lhs == rhs, key
        elif key == "n":
            assert int(lhs) == int(rhs), key
        else:
            assert math.isclose(float(lhs), float(rhs), rel_tol=1e-10, abs_tol=1e-10), (key, lhs, rhs)

def measure_queries(connection, view, aggregates, predicate, expected, row_count):
    result = {}
    for label, suffix in (("full_scan", ""), ("filtered_scan", " WHERE " + predicate)):
        sql = "SELECT " + aggregates + " FROM " + view + suffix

        check_aggregates(expected[label], connection.execute(sql)[0])
        samples = []
        for _ in range(3):
            output, ns = connection.execute(sql)
            check_aggregates(expected[label], output)
            samples.append(ns)
        result[label] = {"samples_ns": samples, "median_ns_per_input_row": statistics.median(samples) / row_count,
                         "output_rows": int(expected[label][0]["n"])}
    return result

def verify(connection, view, names):
    columns = ",".join("cluster_bits(" + name + ") AS " + name for name in names)
    connection.execute("CREATE OR REPLACE TEMP VIEW actual_bits AS SELECT " + columns + " FROM " + view)
    output, _ = connection.execute(
        "SELECT count(*) AS mismatches FROM ((FROM expected_bits EXCEPT ALL FROM actual_bits) "
        "UNION ALL (FROM actual_bits EXCEPT ALL FROM expected_bits))")
    assert int(output[0]["mismatches"]) == 0, "Bitwise complete-tuple multiset mismatch"
    connection.execute("DROP VIEW actual_bits")

def save(report, path):
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", choices=("selected", "full_public"), default="selected")
    parser.add_argument("--datasets", help="Optional comma-separated subset; default all datasets in the scope")
    args = parser.parse_args()
    manifest = WORKSPACE / "data/selected_columns.csv"
    entries = list(csv.DictReader(manifest.open(encoding="utf-8-sig")))
    for entry in entries:
        entry["SourcePath"] = "data/raw/" + entry["Dataset"] + ".csv"
        entry["OutputPath"] = "data/final/" + entry["Dataset"] + ".csv"
    if args.scope in ("selected", "full_public"):
        entries = [e for e in entries if e["Dataset"] in ("Musk", "Crop", "Gas", "CT", "BT")]
    if args.datasets:
        requested = args.datasets.split(",")
        entries = [e for e in entries if e["Dataset"] in requested]
        assert {e["Dataset"] for e in entries} == set(requested), "Unknown dataset"

    path_field = "OutputPath" if args.scope == "selected" else "SourcePath"
    entries.sort(key=lambda e: (WORKSPACE / e[path_field]).stat().st_size)
    check([e['Dataset'] for e in entries], 'final' if args.scope == 'selected' else 'raw')
    out = ROOT / "output" / ("real_" + args.scope + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
    out.mkdir(parents=True)
    report = {"status": "RUNNING", "scope": args.scope, "output_directory": out.relative_to(ROOT).as_posix(),
              "duckdb_version": version().decode(), "duckdb_build_commit": "d591bb1da2",
              "manifest": manifest.relative_to(WORKSPACE).as_posix(), "manifest_sha256": sha256(manifest),
              "parameters": {"page_size": 10000, "k": 100, "seed": 20260929,
                             "compression": "zstd", "compression_level": 3, "threads": 1},
              "protocol": {"input": "all current CSV rows parsed into DOUBLE columns in DuckDB; no resampling",
                           "encoding": "one SQL write per format from the same in-memory typed table; includes file I/O and synchronization; CSV load excluded",
                           "queries": "same-process warm-cache full/filtered aggregates over all columns; one warm-up and three measured calls; C API materialized results",
                           "correctness": "exact binary64 full-tuple multiset with bidirectional EXCEPT ALL for every format; aggregate tolerance only handles summation order",
                           "compression_ratio": "8 * actual_rows * dimension / complete_file_bytes",
                           "limitations": "single write timing, no cold-cache claim, Release -O1, page layouts differ from Parquet; selected scope is not the full original table"},
              "implementation_sha256": {}, "datasets": []}
    for path in [DLL, Path(__file__), ROOT / "scripts/duckdb_local.py",
                 *sorted((ROOT / "src").glob("*.cpp")),
                 *sorted((ROOT / "include").glob("*.hpp")),
                 *sorted((ROOT / "codec").glob("*.cpp")),
                 *sorted((ROOT / "codec").glob("*.h"))]:
        report["implementation_sha256"][str(path.relative_to(ROOT))] = sha256(path)
    print("RESULTS " + str(out), flush=True)
    save(report, out / "report.json")
    for entry in entries:
        name = entry["Dataset"]
        source = WORKSPACE / entry[path_field]
        dim = int(entry["SelectedDim"] if args.scope == "selected" else entry["OriginalDim"])
        names = ["c" + str(c) for c in range(dim)]
        types = "{" + ",".join(literal(c) + ":'DOUBLE'" for c in names) + "}"
        directory = out / name
        directory.mkdir()
        item = {"dataset": name, "input_path": source.relative_to(WORKSPACE).as_posix(), "input_sha256": sha256(source),
                "dimension": dim, "manifest_rows": int(entry["Rows"]),
                "original_column_ids": entry["SelectedColumns"] if args.scope == "selected" else list(range(dim)),
                "formats": {}}
        report["datasets"].append(item)
        print("LOAD " + name, flush=True)
        with Connection() as db:
            _, load_ns = db.execute("CREATE TABLE source AS SELECT * FROM read_csv(" + literal(source) +
                                    ",header=false,auto_detect=false,delim=',',parallel=false,columns=" + types + ")")
            count, _ = db.execute("SELECT count(*) n FROM source")
            n = int(count[0]["n"])
            assert n == int(entry["Rows"]), "Manifest row mismatch"
            item.update(rows=n, typed_bytes=8*n*dim, csv_load_ns=load_ns, manifest_rows_match=n == int(entry["Rows"]))
            bit_columns = ",".join("cluster_bits(" + c + ") AS " + c for c in names)
            db.execute("CREATE TEMP VIEW expected_bits AS SELECT " + bit_columns + " FROM source")
            median, _ = db.execute("SELECT quantile_disc(c0,0.5) threshold FROM source")
            predicate = "c0>" + median[0]["threshold"]
            item["filter"] = predicate
            aggregates = "count(*) n," + ",".join("sum(abs(" + c + ")) s" + str(i) for i, c in enumerate(names))
            expected = {label: db.execute("SELECT " + aggregates + " FROM source" + suffix)[0]
                        for label, suffix in (("full_scan", ""), ("filtered_scan", " WHERE " + predicate))}
            for label, method, compression in (("raw_zstd", "raw", "zstd"), ("acluster", "acluster", "none"),
                                                ("acluster_zstd", "acluster", "zstd"), ("kcluster", "kcluster", "none"),
                                                ("kcluster_zstd", "kcluster", "zstd"), ("parquet_zstd", None, "zstd"),
                                                ("duckdb_native", None, None)):
                print("WRITE " + name + " " + label + " rows=" + str(n) + " dim=" + str(dim), flush=True)
                suffix = ".cluster" if method else ".parquet" if compression else ".duckdb"
                path = directory / (label + suffix)
                if method:
                    options = "FORMAT cluster,METHOD " + literal(method) + ",K 100,PAGE_SIZE 10000,SEED 20260929,COMPRESSION " + literal(compression) + ",COMPRESSION_LEVEL 3"
                    _, write_ns = db.execute("COPY source TO " + literal(path) + " (" + options + ")")
                    view = "read_cluster(" + literal(path) + ")"
                elif compression:
                    _, write_ns = db.execute("COPY source TO " + literal(path) + " (FORMAT parquet,COMPRESSION zstd,COMPRESSION_LEVEL 3,ROW_GROUP_SIZE 10000)")
                    view = "read_parquet(" + literal(path) + ")"
                else:
                    db.execute("ATTACH " + literal(path) + " AS disk")
                    _, write_ns = db.execute("CREATE TABLE disk.data AS SELECT * FROM source; CHECKPOINT disk")
                    db.execute("DETACH disk; ATTACH " + literal(path) + " AS disk (READ_ONLY)")
                    view = "disk.data"
                size = path.stat().st_size
                metrics = {"bytes": size, "compression_ratio": 8*n*dim/size,
                           "write_ns": write_ns, "write_ns_per_row": write_ns/n,
                           "artifact_sha256": sha256(path)}
                if method:
                    records, _ = db.execute("SELECT * FROM cluster_metadata(" + literal(path) + ")")
                    pages = [{k: v if k in ("mode", "fallback_reason", "compression") else int(v)
                              for k, v in record.items()} for record in records]
                    assert sum(p["rows"] for p in pages) == n and len(pages) == (n+9999)//10000
                    metrics.update(pages=pages, encoded_pages=sum(p["mode"] != "raw" for p in pages),
                                   raw_pages=sum(p["mode"] == "raw" for p in pages),
                                   corrections=sum(p["corrections"] for p in pages),
                                   codec_page_bytes=sum(p["codec_page_bytes"] for p in pages),
                                   stored_page_bytes=sum(p["page_bytes"] for p in pages))
                    if label in ("acluster_zstd", "kcluster_zstd"):
                        previous = item["formats"][method]["pages"]
                        fields = ("rows", "mode", "fallback_reason", "references", "corrections", "codec_page_bytes")
                        assert [[p[k] for k in fields] for p in pages] == [[p[k] for k in fields] for p in previous]
                print("VERIFY " + name + " " + label, flush=True)
                verify(db, view, names)
                metrics["bitwise_tuple_multiset"] = "PASS"
                metrics["queries"] = measure_queries(db, view, aggregates, predicate, expected, n)
                if label == "duckdb_native":
                    db.execute("DETACH disk")
                item["formats"][label] = metrics
                save(report, out / "report.json")
                print("PASS " + name + " " + label + " bytes=" + str(size) + " CR=" + format(metrics["compression_ratio"], ".4f") +
                      " write_s=" + format(write_ns/1e9, ".3f"), flush=True)
    report["status"] = "PASS"
    save(report, out / "report.json")
    with (out / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["dataset", "rows", "dimension", "format", "file_bytes", "CR", "write_ns_per_row",
                         "full_scan_median_ns_per_row", "filtered_scan_median_ns_per_row", "raw_pages", "corrections"])
        for dataset in report["datasets"]:
            for label, metrics in dataset["formats"].items():
                writer.writerow([dataset["dataset"], dataset["rows"], dataset["dimension"], label, metrics["bytes"],
                                 metrics["compression_ratio"], metrics["write_ns_per_row"],
                                 metrics["queries"]["full_scan"]["median_ns_per_input_row"],
                                 metrics["queries"]["filtered_scan"]["median_ns_per_input_row"],
                                 metrics.get("raw_pages", ""), metrics.get("corrections", "")])
    print("COMPLETE " + str(out / "report.json"), flush=True)

if __name__ == "__main__":
    main()
