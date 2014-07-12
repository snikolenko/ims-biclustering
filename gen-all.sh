#!/bin/bash

d=$1
echo "Running ims-bicluster for $d..."
./bin/ims-bicluster --mat --input=data/$d --eigens=20
gen-report $d
suffix=$(date '+%y%m%d-%H%M%S')
echo "Making report for $d..."
reportname=$(gen-report $d)
echo "Report written to reports/$reportname"
exit
