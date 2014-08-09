#!/bin/bash

d=$1
mode=$2
echo "Running ims-bicluster for $d..."
./bin/ims-bicluster --mat$mode --input=data/$d --eigens=20
suffix=$(date '+%y%m%d-%H%M%S')
echo "Making report for $d..."
reportname=$(bash gen-report.sh $d $mode)
echo "Report written to reports/$reportname.pdf"
exit
