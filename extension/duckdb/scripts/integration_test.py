import datetime
import json
import os
from pathlib import Path
import statistics
import struct
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / '.build' / ('duckdb.exe' if os.name == 'nt' else 'duckdb')
OUT = ROOT / 'output' / datetime.datetime.now().strftime('integration_%Y%m%d_%H%M%S')
OUT.mkdir(parents=True)
ENV = dict(os.environ)

def literal(path):
    return "'" + str(path).replace("\\", "/").replace("'", "''") + "'"

def run(sql, expect_error=False):
    started = time.perf_counter()
    result = subprocess.run(
        [str(EXE), "-batch", "-json", "-c", "LOAD cluster; SET threads=1; " + sql],
        cwd=ROOT, env=ENV, capture_output=True, text=True, timeout=180,
    )
    duration = time.perf_counter() - started
    if expect_error:
        if result.returncode == 0:
            raise AssertionError("Expected failure: " + sql)
        if "IO Error:" not in result.stderr:
            raise AssertionError("Expected a recoverable I/O error: " + result.stderr)
        return result.stderr.strip(), duration
    if result.returncode:
        raise RuntimeError(result.stderr + "\nSQL: " + sql)
    records = []
    remaining = result.stdout.strip()
    while remaining:
        record, end = json.JSONDecoder().raw_decode(remaining)
        records.append(record)
        remaining = remaining[end:].strip()
    return records, duration

SOURCE = "SELECT (i%17)::DOUBLE*0.25 a,(i%11)::DOUBLE-5 b,(i%7)::DOUBLE*2 c FROM range(20037) t(i)"
report = {"kind": "correctness and smoke comparison; not a paper benchmark", "rows": 20037,
          "columns": 3, "threads": 1, "build": "GCC Release -O1", "files": {}, "checks": []}
expected_query = "SELECT sum(a) sa, sum(b) sb, sum(c) sc FROM source WHERE a>1"
expected, _ = run("CREATE VIEW source AS " + SOURCE + ";" + expected_query)
for method in ("acluster", "kcluster", "raw"):
    path = OUT / (method + ".cluster")
    _, write_seconds = run("COPY (" + SOURCE + ") TO " + literal(path) +
                           " (FORMAT cluster, METHOD '" + method + "', K 16, PAGE_SIZE 10000, SEED 42)")
    view = "read_cluster(" + literal(path) + ")"
    comparison = (
        "CREATE VIEW source AS " + SOURCE + ";"
        "CREATE VIEW lhs AS SELECT cluster_bits(a) a,cluster_bits(b) b,cluster_bits(c) c FROM source;"
        "CREATE VIEW rhs AS SELECT cluster_bits(a) a,cluster_bits(b) b,cluster_bits(c) c FROM " + view + ";"
        "SELECT count(*) mismatches FROM ((FROM lhs EXCEPT ALL FROM rhs) UNION ALL (FROM rhs EXCEPT ALL FROM lhs));"
    )
    rows, _ = run(comparison)
    assert rows[-1][0]["mismatches"] == 0, method
    metadata, _ = run("SELECT * FROM cluster_metadata(" + literal(path) + ")")
    pages = [{key: int(value) if key not in ("mode", "fallback_reason", "compression") else value
              for key, value in page.items()} for page in metadata[-1]]
    assert len(pages) == 3 and sum(p["rows"] for p in pages) == 20037
    if method != "raw":
        assert any(p["mode"] == method for p in pages), "Only RAW pages: " + method
    query, _ = run(expected_query.replace("source", view))
    assert query == expected
    times = [run(expected_query.replace("source", view))[1] for _ in range(3)]
    repeat = OUT / (method + "_repeat.cluster")
    run("COPY (" + SOURCE + ") TO " + literal(repeat) +
        " (FORMAT cluster, METHOD '" + method + "', K 16, PAGE_SIZE 10000, SEED 42)")
    assert path.read_bytes() == repeat.read_bytes(), "Non-deterministic file: " + method
    report["files"][method] = {"bytes": path.stat().st_size, "write_process_seconds": write_seconds,
                                "query_process_median_seconds": statistics.median(times), "pages": pages}
    report["checks"].append(method + ": fresh-process bitwise tuple multiset, aggregates, page totals, deterministic file")

parquet = OUT / "baseline.parquet"
run("COPY (" + SOURCE + ") TO " + literal(parquet) + " (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 10000)")
pq_query = expected_query.replace("source", "read_parquet(" + literal(parquet) + ")")
assert run(pq_query)[0] == expected
report["files"]["parquet_zstd"] = {"bytes": parquet.stat().st_size,
    "query_process_median_seconds": statistics.median(run(pq_query)[1] for _ in range(3))}
database = OUT / "baseline.duckdb"
run("ATTACH " + literal(database) + " AS disk; CREATE TABLE disk.source AS " + SOURCE + "; CHECKPOINT disk;")
db_query = "ATTACH " + literal(database) + " AS disk (READ_ONLY);" + expected_query.replace("source", "disk.source")
assert run(db_query)[0] == expected
report["files"]["duckdb_native"] = {"bytes": database.stat().st_size,
    "query_process_median_seconds": statistics.median(run(db_query)[1] for _ in range(3))}

original = (OUT / "acluster.cluster").read_bytes()
header_length = struct.unpack_from("<Q", original, 8)[0]
for label, contents in (
    ("truncated", original[:-1]), ("header", original[:16] + bytes([original[16] ^ 1]) + original[17:]),
    ("page", original[:16 + header_length + 30] + bytes([original[16 + header_length + 30] ^ 1]) + original[17 + header_length + 30:]),
    ("trailer", original[:-1] + bytes([original[-1] ^ 1])),
):
    path = OUT / ("corrupt_" + label + ".cluster")
    path.write_bytes(contents)

    message, _ = run("CREATE TEMP TABLE should_fail AS SELECT * FROM read_cluster(" + literal(path) + ")",
                     expect_error=True)
    report["checks"].append(label + " rejected")

empty = OUT / "empty.cluster"
run("COPY (SELECT 1::DOUBLE a WHERE false) TO " + literal(empty) + " (FORMAT cluster)")
assert run("SELECT count(*) n FROM read_cluster(" + literal(empty) + ")")[0][-1][0]["n"] == 0
report["checks"].append("Empty file schema and zero rows")
report["limitations"] = ["DOUBLE only; immutable local files; sequential scan; no predicate/projection pushdown",
    "Process timings include CLI startup; tiny synthetic comparison is not evidence of general performance",
    "Parquet row groups and cluster pages are different units; DuckDB file size includes database allocation overhead",
    "Corruption tests force complete result consumption with CTAS"]
(OUT / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps({"status": "PASS", "report": str(OUT / "report.json"),
                  "bytes": {k: v["bytes"] for k, v in report["files"].items()}}, ensure_ascii=False))
