#!/bin/bash
for kmean in {18..18}
do
for seed in {1..100}
do
echo working on $kmean\, seed = $seed
python script/Analyze4.py $kmean $seed >> test.txt
done
done
