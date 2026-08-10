import argparse
import csv
import datetime
from pathlib import Path
import subprocess
from check_data import ROOT, check
from summarize import summarize

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('experiment', choices=['end-to-end', 'optimal', 'order'])
    parser.add_argument('--datasets', help='Comma-separated dataset names')
    parser.add_argument('--method', default='all')
    parser.add_argument('--runs', type=int)
    parser.add_argument('--warmups', type=int, default=3)
    parser.add_argument('--layout', choices=['original', 'page-sort'], default='original')
    parser.add_argument('--output', type=Path, help='Relative to encoding_benchmark')
    args = parser.parse_args()
    defaults = {'end-to-end':['Musk','Crop','Gas','CT','BT'], 'optimal':['SSD','BC','PI'], 'order':['BC','BT']}
    names = args.datasets.split(',') if args.datasets else defaults[args.experiment]
    manifest = list(csv.DictReader((ROOT / 'data/selected_columns.csv').open()))
    valid = {row['Dataset'] for row in manifest}
    if not names or not set(names) <= valid:
        parser.error('Unknown dataset')
    if args.experiment == 'end-to-end' and not set(names) <= set(defaults['end-to-end']):
        parser.error('End-to-end requires a multidimensional public dataset')
    if args.experiment == 'optimal' and not set(names) <= set(defaults['optimal']):
        parser.error('Optimal supports SSD, BC and PI')
    runs = args.runs if args.runs is not None else (3 if args.experiment == 'order' else 10)
    if runs < 1 or args.warmups < 0:
        parser.error('Invalid repeat count')
    out = args.output or Path('output') / (args.experiment + '_' + datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
    if out.is_absolute() or '..' in out.parts:
        parser.error('Use a repository-relative output directory')
    out = ROOT / out
    if out.exists():
        parser.error('Use a fresh output directory')
    scope = 'raw' if args.experiment == 'end-to-end' else 'final'
    check(names, scope)
    binary = 'order_benchmark' if args.experiment == 'order' else 'table_benchmark'
    executable = next((ROOT / 'build' / (binary+suffix) for suffix in ['','.exe'] if (ROOT / 'build' / (binary+suffix)).is_file()), None)
    if not executable:
        parser.error('Build the benchmark first')
    for name in names:
        if args.experiment == 'order':
            command = [str(executable), 'data/final/'+name+'.csv', str(out/name), args.layout, 'both', str(runs)]
        else:
            command = [str(executable), '--experiment', args.experiment, '--dataset', name,
                       '--output', str(out/name if args.experiment=='optimal' else out), '--runs', str(runs),
                       '--warmups', str(args.warmups), '--method', args.method]
        subprocess.run(command, cwd=ROOT, check=True)
    summarize(out, args.experiment)
    print('Experiment complete')

if __name__ == '__main__':
    main()
