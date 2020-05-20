from collections import defaultdict
import sys
from glob import glob

examplestats = defaultdict(list)
filestats = {}

fns_glob = glob('*csv')
for fn in fns_glob:
    csv = open(fn)
    next(csv)
    total_time = 0
    total_solved = 0
    for line in csv:
        cols = line.split(',')
        if len(cols) == 4:
            example = cols[0]
            td = int(cols[1])
            time = float(cols[2])
            if (td > -1):
                total_time += time
                total_solved += 1
                examplestats[example].append((fn, td, time))
    filestats[fn] = (total_time, total_solved)

print('Example,Treedepth,{}'.format(','.join([fn.split('/')[-1].replace("computed_treedepths_","").replace(".csv","") for fn in fns_glob])), file=sys.stderr)
for example in sorted(examplestats.keys()):
    stats = examplestats[example]
    td_min = 99999
    time_min = 9999
    fns = {}
    for fn, td, time in stats:
        td_min = min(td, td_min)
        fns[fn] = time
    for fn, td, time in stats:
        if td != td_min:
            print('Error!? In {} the td was {} whereas the minimal td was {}.'.
                  format(fn, td, td_min))

    times = []
    for fn in fns_glob:
        if fn in fns:
            times.append(fns[fn])
        else:
            times.append(-1)
    #fns.sort(key=lambda fn: -filestats[fn][1])
    print('{},{:<3},{}'.format(example, td_min, ','.join([str(time).rjust(6) for time in times])), file=sys.stderr)

print('')
print('File\t\t\t\t\tTotal Solved\tAverage time')
fns_glob.sort(key=lambda fn: -filestats[fn][1])
for fn in fns_glob:
    total_time, total_solved = filestats[fn]
    if total_solved:
        print('{:<30}\t\t{}\t\t{}'.format(fn.split('/')[-1].replace("computed_treedepths_","").replace(".csv",""), total_solved,
                                      total_time / total_solved))
print('\nTotal solved:', len(examplestats))
