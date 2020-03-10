from glob import glob

fns = glob('*csv')
for fn in fns:
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
    if total_solved:
        print('File {} solved a total of {} cases with average solve time of {}'.format(fn, total_solved, total_time / total_solved))

