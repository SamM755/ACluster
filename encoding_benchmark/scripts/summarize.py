import argparse
import csv
from collections import defaultdict
from pathlib import Path
import statistics

def summarize(root, experiment):
    groups = defaultdict(list)
    pattern = '*_runs.csv' if experiment == 'end-to-end' else 'runs.csv'
    if experiment == 'order':
        rows = []
        for path in sorted(root.glob('*/summary.csv')):
            rows.extend(csv.DictReader(path.open()))
        if not rows:
            raise RuntimeError('No order results')
        columns = list(rows[0])
    else:
        for path in sorted(root.rglob(pattern)):
            if 'warmup' in path.name:
                continue
            for row in csv.DictReader(path.open()):
                if row['Verified'].lower() not in ['true','1']:
                    raise RuntimeError('Unverified measurement')
                groups[(row['Dataset'],row['Method'])].append(row)
        if not groups:
            raise RuntimeError('No measurements')
        rows=[]
        for (dataset,method), values in groups.items():
            byte_key='TotalBytes' if experiment=='end-to-end' else 'Bytes'
            if len({row[byte_key] for row in values})!=1:
                raise RuntimeError('Non-deterministic stream size')
            result={'Dataset':dataset,'Method':method,'Runs':len(values),'Bytes':int(values[0][byte_key]),'CR':float(values[0]['CR'])}
            for metric in ['EncodeNsPerRow','DecodeNsPerRow']:
                samples=[float(row[metric]) for row in values]
                result[metric]=statistics.mean(samples)
                result[metric+'Std']=statistics.stdev(samples) if len(samples)>1 else 0
            result['Verified']='true'
            rows.append(result)
        columns=list(rows[0])
    with (root/'summary.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns)
        writer.writeheader();writer.writerows(rows)

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('experiment',choices=['end-to-end','optimal','order'])
    parser.add_argument('directory',type=Path)
    args=parser.parse_args()
    summarize(args.directory,args.experiment)
