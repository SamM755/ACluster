import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def check(names=None, scope=None):
    hashes = json.loads((ROOT / 'data/checksums.json').read_text())
    for relative, expected in hashes.items():
        path = ROOT / relative
        if names and path.stem not in names:
            continue
        if scope and path.parent.name != scope:
            continue
        if not path.is_file():
            raise RuntimeError('Missing dataset: ' + relative)
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            chunk = stream.read(1024 * 1024)
            if chunk.startswith(b'version https://git-lfs.github.com/spec/v1'):
                raise RuntimeError('Git LFS pointer: run git lfs pull for ' + relative)
            while chunk:
                digest.update(chunk)
                chunk = stream.read(1024 * 1024)
        if digest.hexdigest() != expected:
            raise RuntimeError('Dataset checksum mismatch: ' + relative)

if __name__ == '__main__':
    check()
    print('All dataset checksums passed')
