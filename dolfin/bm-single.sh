#!/bin/bash


for i in `seq 1 5`; do
{
  echo Run $i;
  python advdiff.py $1 > stats/output$1-run$i.txt;
}
done;
