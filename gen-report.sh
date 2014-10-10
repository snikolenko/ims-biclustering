#!/bin/bash

d=$1
mattype=$2
suffix=$(date '+%y%m%d-%H%M%S')
if [ -z "$3" ]
  then
    M=10
  else
  	M=$3
fi
if [ -z "$4" ]
  then
    N=10
  else
  	N=$4
fi
python python/visualize.py -x data/$d -y $mattype -e data/$d.val.csv -v data/$d.vec.csv -c data/$d.coords.csv -N $N -M $M -t $d
cd reports/latex && pdflatex tmp >/dev/null && pdflatex tmp >/dev/null && mv tmp.pdf ../$d-$suffix.pdf && cd ../..
echo $d-$suffix
exit
