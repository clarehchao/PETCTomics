#!/bin/bash

for n in {1..10}
do
   num=$RANDOM
   a=$(echo "scale=5; $num / 32767" | bc)
   b=$(echo "scale=5; $a * 75" | bc)
   c=$(echo "scale=5; $b + 280" | bc)
   logfname="ELNrunlog$n.txt"
   ./Radiomics/classification_cv.py param_file/ELN_TGBinary.json $n > $logfname &
   echo "submitted Run#$n!"
   echo "sleep for $c seconds!"
   sleep $c
done

echo "All jobs submitted!"

