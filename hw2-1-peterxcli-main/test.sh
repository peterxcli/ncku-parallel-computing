#!/bin/bash

test=$1

make
mpiexec -f hosts -n 2 ./hw2_1 < input/filename/"$test".txt > tmp.out
	diff tmp.out answer/"$test".out
	rm tmp.out
