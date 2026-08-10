import argparse
import ctypes
import json
from pathlib import Path
import struct

from benchmark_real import check_aggregates, verify
from duckdb_local import Connection, literal, ROOT
from validate_zstd_files import seal

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    report = json.loads((args.results / "report.json").read_text(encoding="utf-8"))
    assert report["status"] == "PASS"
    scratch = args.results / "format_validation"
    scratch.mkdir()
    checks = []
    with Connection() as db:
        for item in report["datasets"]:
            print("REOPEN " + item["dataset"], flush=True)
            names = [f"c{c}" for c in range(item["dimension"])]
            expected_path = ROOT / item["formats"]["raw_zstd"]["runs"][0]["path"]
            columns = ",".join("cluster_bits(" + name + ") AS " + name for name in names)
            db.execute("CREATE OR REPLACE TEMP VIEW expected_bits AS SELECT " + columns +
                       " FROM read_cluster(" + literal(expected_path) + ")")
            for label, metrics in item["formats"].items():
                path = ROOT / metrics["runs"][0]["path"]
                view = ("read_parquet(" if label == "parquet_zstd" else "read_cluster(") + literal(path) + ")"
                check_aggregates(item["expected_aggregates"]["full_scan"],
                                 db.execute("SELECT " + item["aggregates"] + " FROM " + view)[0])
                verify(db, view, names)
                checks.append(item["dataset"] + " " + label + ": fresh-process complete-tuple check PASS")
        original_path = scratch / "hybrid_uncompressed.cluster"
        db.execute("COPY (SELECT ((i*37)%7)::DOUBLE a, (i%997)::DOUBLE/100 b, ((i*37)%7)::DOUBLE*2 c, "
                   "(i%13)::DOUBLE d FROM range(1031) t(i)) TO " + literal(original_path) +
                   " (FORMAT cluster,METHOD acluster_alp,CLUSTER_COLUMNS '2,0',K 8,COMPRESSION none)")
        original = original_path.read_bytes()
        footer = struct.unpack_from("<Q", original, len(original)-24)[0]
        offset, length = struct.unpack_from("<QQ", original, footer+16)
        page = original[offset:offset+length]
        assert struct.unpack_from("<I", page, 4)[0] == 3
        selected = struct.unpack_from("<I", page, 28)[0]
        assert selected == 2
        cluster_len_position = 32 + selected * 4
        cluster_length = struct.unpack_from("<Q", page, cluster_len_position)[0]
        column_position = cluster_len_position + 8 + cluster_length
        assert page[column_position] == 1
        alp_position = column_position + 9
        variants = []
        invalid_column = bytearray(page[:-8])
        struct.pack_into("<I", invalid_column, 32+4, struct.unpack_from("<I", page, 32)[0])
        variants.append(("duplicate_selected_column", seal(invalid_column)))
        invalid_alp_count = bytearray(page[:-8])
        struct.pack_into("<I", invalid_alp_count, alp_position, 1032)
        variants.append(("alp_total_count", seal(invalid_alp_count)))
        invalid_block = bytearray(page[:-8])
        struct.pack_into("<H", invalid_block, alp_position+5, 0)
        variants.append(("alp_zero_block_count", seal(invalid_block)))
        damaged = bytearray(page)
        damaged[45] ^= 1
        variants.append(("hybrid_page_checksum", damaged))
        for label, replacement in variants:
            path = scratch / (label + ".cluster")
            path.write_bytes(original[:offset] + replacement + original[offset+length:])
            try:
                db.execute("SELECT count(*) FROM read_cluster(" + literal(path) + ")")
            except RuntimeError as error:
                checks.append(label + ": rejected")
            else:
                raise AssertionError("Corrupt hybrid file accepted: " + label)
            assert db.execute("SELECT 42 ok")[0][0]["ok"] == "42"
    result = {"status": "PASS", "checks": checks}
    (scratch / "validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
