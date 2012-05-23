#!/bin/bash

rm -rf stats
mkdir stats

for i in `seq  1 12`; do
{
  echo Running $i
  ./bm-single.sh $i
}
done;
