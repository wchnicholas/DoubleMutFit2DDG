#!/bin/bash
for n in {1..100}
do
grep -P "^$n\ " $1 | sort -k 3 -n | head -1
done
