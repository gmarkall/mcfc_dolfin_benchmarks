#!/bin/bash

for i in `seq  1 12`; do
{
  echo Running $i
  ./profile.sh $i
}
done;
