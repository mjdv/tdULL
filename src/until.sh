#!/bin/bash

INPUTDIR="../input/exact/"
OUTPUTDIR="../output/"
for f in ${INPUTDIR}*; do
    graph="${f##*/exact_}"
    graph="${graph%.*}"

    ./main ${INPUTDIR}exact_${graph}.gr ${OUTPUTDIR}exact_${graph}.tree;

    if [ "${graph}" = "$1" ] ; then
      break;
    fi
done  


