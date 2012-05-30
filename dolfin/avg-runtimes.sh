#!/bin/bash

rm -f stats/times.txt

echo -e "Cores\tTime" > stats/times.txt

for i in `seq 1 12`; do
{
  a=`cat stats/output$i-run*.txt | grep Timestepping | awk 'BEGIN { s=0 } {s=s+$5} END {print s/NR}'`;
  echo -e "$i\t$a" >> stats/times.txt;
}
done;
