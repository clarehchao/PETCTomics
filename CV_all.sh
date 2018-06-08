#!/bin/bash

num=$RANDOM
a=$(echo "scale=5; $num / 32767" | bc)
b=$(echo "scale=5; $a * 75" | bc)
c=$(echo "scale=5; $b + 280" | bc)
./Radiomics/classification_cv.py param_file/L1LiblinearLogReg_TGBinary.json > l1liblinearRunlog.txt &
echo "sleep for $c seconds!"
sleep $c

num=$RANDOM
a=$(echo "scale=5; $num / 32767" | bc)
b=$(echo "scale=5; $a * 75" | bc)
c=$(echo "scale=5; $b + 280" | bc)
./Radiomics/classification_cv.py param_file/L1SagaLogReg_TGBinary.json > l1SagaRunlog.txt &
echo "sleep for $c seconds!"
sleep $c

num=$RANDOM
a=$(echo "scale=5; $num / 32767" | bc)
b=$(echo "scale=5; $a * 75" | bc)
c=$(echo "scale=5; $b + 280" | bc)
./Radiomics/classification_cv.py param_file/L2LiblinearLogReg_TGBinary.json > l2LiblinearRunlog.txt &
echo "sleep for $c seconds!"
sleep $c

num=$RANDOM
a=$(echo "scale=5; $num / 32767" | bc)
b=$(echo "scale=5; $a * 75" | bc)
c=$(echo "scale=5; $b + 280" | bc)
./Radiomics/classification_cv.py param_file/L2Newtoncg_TGBinary.json > l2NewtonRunlog.txt &
echo "sleep for $c seconds!"
sleep $c

num=$RANDOM
a=$(echo "scale=5; $num / 32767" | bc)
b=$(echo "scale=5; $a * 75" | bc)
c=$(echo "scale=5; $b + 280" | bc)
./Radiomics/classification_cv.py param_file/L2LbfgsLogReg_TGBinary.json > l2lbfgsRunlog.txt &
echo "sleep for $c seconds!"
sleep $c

num=$RANDOM
a=$(echo "scale=5; $num / 32767" | bc)
b=$(echo "scale=5; $a * 75" | bc)
c=$(echo "scale=5; $b + 280" | bc)
./Radiomics/classification_cv.py param_file/L2SagLogReg_TGBinary.json > l2SagRunlog.txt &
echo "sleep for $c seconds!"
sleep $c

num=$RANDOM
a=$(echo "scale=5; $num / 32767" | bc)
b=$(echo "scale=5; $a * 75" | bc)
c=$(echo "scale=5; $b + 280" | bc)
./Radiomics/classification_cv.py param_file/SVM_TGBinary.json > SVMRunlog.txt &
echo "sleep for $c seconds!"
sleep $c

num=$RANDOM
a=$(echo "scale=5; $num / 32767" | bc)
b=$(echo "scale=5; $a * 75" | bc)
c=$(echo "scale=5; $b + 280" | bc)
./Radiomics/classification_cv.py param_file/RandomForest_TGBinary.json > RFRunlog.txt &
echo "sleep for $c seconds!"
sleep $c


















