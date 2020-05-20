#!/bin/bash
INPUTDIR="../input/exact/"
OUTPUTDIR="../output/exact/"
for N in $(seq 1 28);
do
    rm main.o
    make BATCHMETHOD=B$N > /dev/null
    CSVFILE="timings/B$N.csv"
    cat /dev/null > ${CSVFILE}

    echo "Dumping timings in $CSVFILE";
    for f in exact_001.gr exact_003.gr exact_005.gr exact_007.gr exact_009.gr exact_011.gr exact_013.gr exact_015.gr exact_017.gr exact_019.gr exact_021.gr exact_023.gr exact_025.gr exact_027.gr exact_029.gr exact_031.gr exact_033.gr exact_035.gr exact_037.gr exact_039.gr exact_041.gr exact_043.gr exact_045.gr exact_047.gr exact_049.gr exact_051.gr exact_053.gr exact_055.gr exact_057.gr exact_059.gr exact_061.gr exact_063.gr exact_065.gr exact_067.gr exact_069.gr exact_071.gr exact_073.gr exact_077.gr exact_079.gr exact_081.gr exact_083.gr exact_085.gr exact_087.gr exact_089.gr exact_091.gr exact_093.gr exact_095.gr exact_097.gr exact_099.gr exact_103.gr exact_105.gr exact_107.gr exact_109.gr exact_111.gr exact_113.gr exact_115.gr exact_117.gr exact_123.gr exact_125.gr exact_127.gr exact_133.gr exact_137.gr exact_141.gr exact_143.gr exact_145.gr exact_151.gr exact_153.gr exact_157.gr exact_159.gr exact_161.gr exact_165.gr exact_173.gr exact_177.gr exact_181.gr exact_185.gr exact_189.gr
    do
        printf "\nGraph ${f}\n"
        timeout 5m ./main < ${INPUTDIR}${f} > /dev/null 2> /tmp/output$$ && printf "${f}," >> $CSVFILE && tail -1 /tmp/output$$ >> $CSVFILE;
        sed -e '$d' /tmp/output$$
    done
done
