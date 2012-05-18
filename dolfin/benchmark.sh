#!/bin/bash

rm -rf output
mkdir output

for i in `seq 1 6`; do
{
  echo "Running with $i procs...";
  mpirun -n $i python run.py | tee output/$i.txt
}
done;
