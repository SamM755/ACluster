import argparse
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', help='Optional local DuckDB Git repository')
    args = parser.parse_args()
    dependency = json.loads((ROOT / 'dependencies.json').read_text())['duckdb']
    target = ROOT / '.deps/duckdb'
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        if args.source:
            subprocess.run(['git', 'clone', '--no-hardlinks', '--no-checkout', str(Path(args.source).resolve()), str(target)], check=True)
        else:
            subprocess.run(['git', 'init', str(target)], check=True)
            subprocess.run(['git', '-C', str(target), 'remote', 'add', 'origin', dependency['url']], check=True)
            subprocess.run(['git', '-C', str(target), 'fetch', '--depth', '1', 'origin', dependency['commit']], check=True)
        subprocess.run(['git', '-C', str(target), 'checkout', '--detach', dependency['commit']], check=True)
    commit = subprocess.check_output(['git', '-C', str(target), 'rev-parse', 'HEAD'], text=True).strip()
    if commit != dependency['commit']:
        raise RuntimeError('Unexpected DuckDB commit')
    extension = target / 'extension/cluster'
    extension.mkdir(parents=True, exist_ok=True)
    for folder in ['src', 'include', 'codec']:
        shutil.copytree(ROOT / folder, extension / folder, dirs_exist_ok=True)
    shutil.copy2(ROOT / 'CMakeLists.txt', extension / 'CMakeLists.txt')
    tests = target / 'test/sql/cluster'
    tests.mkdir(parents=True, exist_ok=True)
    for path in (ROOT / 'tests').glob('*.test'):
        shutil.copy2(path, tests / path.name)
    print('DuckDB source ready: ' + commit)

if __name__ == '__main__':
    main()
