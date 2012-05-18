#!/bin/bash

rm -rf stats$1
mkdir stats$1

mpirun -n $1 likwid-perfctr -c L:0 -g MEM -o stats$1/advdiff-mem_%p.txt python advdiff.py > stats$1/output1.txt
mpirun -n $1 likwid-perfctr -c L:0 -g FLOPS_DP -o stats$1/advdiff-flops_%p.txt python advdiff.py > stats$1/output2.txt
