#!/bin/bash

#This script is for check bed files, including:
for i in *bed
do
	#1. remove scaffolds and sort 
	grep -vE "GL|JH|KI|chrPt|KN|KZ|Scaffold|chr2110|chrrDNA|ctg" ${i} | sort --parallel=20 -k1,1 -k2,2n > ${i}_tmp && mv ${i}_tmp ${i}
	#2. malformated bed lines $2>$3 | $2 < 0 | $3 < 0
	nohup python ~/Zhang_Scripts/Zhang/check_bedfile_valid.py -b ${i} &
done