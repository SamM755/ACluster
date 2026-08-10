import argparse
import json
from pathlib import Path
import struct

from duckdb_local import Connection, literal

def seal(body):
    value = 14695981039346656037
    for byte in body:
        value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return body + struct.pack("<Q", value)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    report = json.loads((args.results / "report.json").read_text(encoding="utf-8"))
    dataset = report["datasets"][0]
    directory = args.results / dataset["dataset"]
    scratch = args.results / "format_validation"
    scratch.mkdir()
    checks = []
    with Connection() as db:
        for method in ("raw_zstd", "acluster", "acluster_zstd", "kcluster", "kcluster_zstd"):
            path = directory / (method + ".cluster")
            result, _ = db.execute("SELECT count(*) n FROM read_cluster(" + literal(path) + ")")
            assert int(result[0]["n"]) == dataset["rows"]
            checks.append(method + ": reopened in a new process")
        original = (directory / "acluster_zstd.cluster").read_bytes()
        assert original[:8] == b"CLSTR002"
        header_length = struct.unpack_from("<Q", original, 8)[0]
        footer_offset = struct.unpack_from("<Q", original, len(original) - 24)[0]
        page_offset, page_length = struct.unpack_from("<QQ", original, footer_offset + 16)
        assert page_offset == 16 + header_length
        page = bytearray(original[page_offset:page_offset+page_length])
        altered_length = bytearray(page[:-8])
        struct.pack_into("<Q", altered_length, 0, struct.unpack_from("<Q", page)[0]+1)
        invalid_frame = bytearray(page[:-8])
        invalid_frame[8] ^= 255
        def replaced_page(body):
            return original[:page_offset] + seal(body) + original[page_offset+page_length:]
        damaged_page = bytearray(original)
        damaged_page[page_offset+12] ^= 1
        for label, data in (
            ("truncated", original[:-1]),
            ("header_checksum", original[:16] + bytes([original[16] ^ 1]) + original[17:]),
            ("page_checksum", damaged_page),
            ("decoded_length", replaced_page(altered_length)),
            ("invalid_zstd_frame", replaced_page(invalid_frame)),
            ("footer_checksum", original[:footer_offset] + bytes([original[footer_offset] ^ 1]) + original[footer_offset+1:]),
        ):
            path = scratch / (label + ".cluster")
            path.write_bytes(data)
            try:
                db.execute("SELECT count(*) n FROM read_cluster(" + literal(path) + ")")
            except RuntimeError as error:
                assert "IO Error:" in str(error), str(error)
                checks.append(label + ": rejected")
            else:
                raise AssertionError("Corrupt file accepted: " + label)
            assert db.execute("SELECT 42 ok")[0][0]["ok"] == "42"
    result = {"status": "PASS", "checks": checks}
    (scratch / "validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result))

if __name__ == "__main__":
    main()
