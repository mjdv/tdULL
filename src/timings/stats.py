from collections import defaultdict
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
            time = int(cols[2])
            if (td > -1):
                total_time += time
                total_solved += 1
                examplestats[example].append((fn, td, time))
    filestats[fn] = (total_time, total_solved)

print('Example\tTreedepth\tSolved by')
for example in sorted(examplestats.keys()):
    stats = examplestats[example]
    td_min = 99999
    time_min = 9999
    fns = []
    for fn, td, time in stats:
        td_min = min(td, td_min)
        fns.append(fn)
    for fn, td, time in stats:
        if td != td_min:
            print('Error!? In {} the td was {} whereas the minimal td was {}.'.
                  format(fn, td, td_min))

    fns.sort(key=lambda fn: -filestats[fn][1])
    print('{}\t{}\t{}'.format(example, td_min, fns))

print('')
print('File\t\t\t\t\tTotal Solved\tAverage time')
fns_glob.sort(key=lambda fn: -filestats[fn][1])
for fn in fns_glob:
    total_time, total_solved = filestats[fn]
    if total_solved:
        print('{}\t\t{}\t\t{}'.format(fn, total_solved,
                                      total_time / total_solved))
