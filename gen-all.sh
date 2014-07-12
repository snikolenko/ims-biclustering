#!/bin/bash

d=$1
echo "Running ims-bicluster for $d..."
./bin/ims-bicluster --mat --input=data/$d --eigens=20
suffix=$(date '+%y%m%d-%H%M%S')
echo "Making report for $d..."
reportname=$(bash gen-report.sh $d)
echo "Report written to reports/$reportname.pdf"
exit
