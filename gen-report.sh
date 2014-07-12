#!/bin/bash

d=$1
suffix=$(date '+%y%m%d-%H%M%S')
python python/visualize.py -e data/$d.val.csv -v data/$d.vec.csv -c data/$d.coords.csv -N 10 -M 10 -t $d
cd reports/latex && pdflatex tmp >/dev/null && pdflatex tmp >/dev/null && mv tmp.pdf ../$d-$suffix.pdf && cd ../..
echo $d-$suffix
exit
