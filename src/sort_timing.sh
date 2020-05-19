#!/bin/bash
for N in $(seq 1 15);
do
    rm treedepth_test.o
    make SORTMETHOD=SORT$N > /dev/null
    echo "SORT$N" > timings/SORT$N.csv
    ./treedepth_test > timings/SORT$N.csv 2> /dev/null
    printf "SORT$N total time took " && tail -1 timings/SORT$N.csv
done
