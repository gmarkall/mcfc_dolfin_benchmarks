#!/bin/bash

rm -rf stats
mkdir stats

echo -e "Cores\tMemBW\tRemBW\tFLOPS" > stats/average.csv
echo -e "Cores\tMemBW\tRemBW\tFLOPS" > stats/total.csv
echo -e "Cores\tTime" > stats/time.csv

for i in `seq 1 12`; do 
{
  avgmembw=`cat stats$i/membw.txt | awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}'`;
  avgrembw=`cat stats$i/rembw.txt | awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}'`;
  avgflops=`cat stats$i/flops.txt | awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}'`;
  echo -e "$i\t${avgmembw}\t${avgrembw}\t${avgflops}" >> stats/average.csv;

  membw=`cat stats$i/membw.txt | awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}'`;
  rembw=`cat stats$i/rembw.txt | awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}'`;
  flops=`cat stats$i/flops.txt | awk 'BEGIN{s=0;}{s=s+$1;}END{print s;}'`;
  echo -e "$i\t${membw}\t${rembw}\t${flops}" >> stats/total.csv;

  t=`cat stats$i/output1.txt | grep Timestepping | awk '{print $5}'`;
  echo -e "$i\t$t" >> stats/time.csv;
}
done;
