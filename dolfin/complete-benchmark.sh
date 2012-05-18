#!/bin/bash

echo Running all
./runall.sh
echo summarising all
./sumall.sh
echo joining all
./join.sh

echo plotting
gnuplot average.gplt
gnuplot total.gplt
gnuplot runtime.gplt
