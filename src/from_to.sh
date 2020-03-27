#!/bin/bash
FROM=001
UNTIL=199

if [[ "$2" != "" ]]; then
  FROM=$1
  UNTIL=$2
else
  if [[ "$1" != "" ]]; then
    UNTIL=$1
  fi
fi

INPUTDIR="../input/exact/"
OUTPUTDIR="../output/exact/"
GO=0


echo "Computing treedepth of graphs $FROM until $UNTIL"

for f in ${INPUTDIR}*; do
    graph="${f##*/exact_}"
    graph="${graph%.*}"
    if [ "${graph}" = $FROM ] ; then
      GO=1
    fi

    if [ "$GO" == 1 ] ; then
      printf "\nGraph ${graph}\n"
      timeout 11m ./main ${INPUTDIR}exact_${graph}.gr ${OUTPUTDIR}exact_${graph}.tree &&
      ./verify ${INPUTDIR}exact_${graph}.gr ${OUTPUTDIR}exact_${graph}.tree;
    fi

    if [ "${graph}" = "$UNTIL" ] ; then
      break;
    fi
done
