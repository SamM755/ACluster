#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def parse_columns(value: str) -> list[int]:
    return [int(item) for item in value.split(";") if item]


def project(source: Path, destination: Path, columns: list[int]) -> int:
    destination.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with source.open("r", newline="") as source_file, destination.open("w", newline="") as output_file:
        reader = csv.reader(source_file)
        writer = csv.writer(output_file, lineterminator="\n")
        for row in reader:
            if not row or not any(value.strip() for value in row):
                continue
            if max(columns) >= len(row):
                raise RuntimeError(f"{source}: row {rows + 1} has only {len(row)} columns")
            writer.writerow([row[column] for column in columns])
            rows += 1
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Project the selected benchmark columns.")
    parser.add_argument("--data-root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--manifest", type=Path, default=None)
    args = parser.parse_args()

    root = args.data_root.resolve()
    manifest = args.manifest.resolve() if args.manifest else root / "selected_columns.csv"
    with manifest.open("r", newline="") as manifest_file:
        for entry in csv.DictReader(manifest_file):
            name = entry["Dataset"]
            expected_rows = int(entry["Rows"])
            columns = parse_columns(entry["SelectedColumns"])
            written = project(root / "raw" / f"{name}.csv", root / "final" / f"{name}.csv", columns)
            if written != expected_rows:
                raise RuntimeError(f"{name}: wrote {written} rows, expected {expected_rows}")
            print(f"{name}: {written} rows, {len(columns)} columns")


if __name__ == "__main__":
    main()
