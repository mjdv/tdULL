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
	  echo "Calculating treedepth for ${INPUTDIR}exact_${graph}.tree" > /tmp/output
      timeout 31m ./main < ${INPUTDIR}exact_${graph}.gr > ${OUTPUTDIR}exact_${graph}.tree 2>> /tmp/output &&
      printf "exact_${graph}.gr," 1>&2 && tail -1  /tmp/output 1>&2 &&
      #./verify ${INPUTDIR}exact_${graph}.gr ${OUTPUTDIR}exact_${graph}.tree;
      sed -e '$d' /tmp/output
    fi

    if [ "${graph}" = "$UNTIL" ] ; then
      break;
    fi
done
