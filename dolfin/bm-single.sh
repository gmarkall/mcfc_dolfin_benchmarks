#!/bin/bash


mpirun -n $1 --bind-to-socket python advdiff.py > stats/output$1.txt;
