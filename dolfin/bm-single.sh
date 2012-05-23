#!/bin/bash


for i in `seq 1 5`; do
{
  echo Run $i;
  mpirun -n $1 --bind-to-socket python advdiff.py > stats/output$1-run$i.txt;
}
done;
