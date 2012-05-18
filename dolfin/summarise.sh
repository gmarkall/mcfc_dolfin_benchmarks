#!/bin/bash

cat stats$1/advdiff-mem_*.txt | grep "Memory bandwidth" | awk '{ print $6 }' > stats$1/membw.txt
cat stats$1/advdiff-mem_*.txt | grep "Remote BW" | awk '{ print $6 }' > stats$1/rembw.txt
cat stats$1/advdiff-flops_*.txt | grep "DP MFlops" | awk '{ print $5 }' > stats$1/flops.txt

