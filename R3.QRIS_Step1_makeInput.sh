#!/bin/bash
##########################################################################
#This script is for preparing input data for GLM 
#Issue report on Hughiez047@gmail.com
#Copyright (c) 2020 __YenLab@SKLEH__. All rights reserved.
##########################################################################
source /public/home/zhy/scripts/0.utilities.sh

for i in *bed
do
	base=`basename ${i} .bed`
	
	#Get species information for current sample
	id=`echo ${base} | sed "s/_.*//"`
	Configuration_info ${id} 
	echo "--->>> Start Process ${i} for ${GENOME_NAME} ..."

	sample_size=10000
	Prefix=${base}_shuf${sample_size}

	#Skip current sample if has been processed
	if [[ -f ${ft_Prefix}_ext25bp.bed ]]; then
		continue
	fi

	bedtools random -l 51 -n 1000000 -g ${GENOME_SIZE_chrM} > ${base}_random.bed

	#Handle on negative sites file
	bedtools getfasta -fi ${FASTA_FILE} -bed ${base}_random.bed \
	| awk -v RS=">" '!/N/{printf $0RT}' | grep ">" | sed 's/[:-]/\t/g' | sed 's/>//' \
	| awk -v FS="\t" -v OFS="\t" '{if(($3-$2) == 51) {print $0}}' \
	| awk '!a[$1$2$3]++' | head -${sample_size} | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,"0","."}' > ${sample_size}_uncut_ext25bp.bed
	bedtools getfasta -fi ${FASTA_FILE} -bed ${sample_size}_uncut_ext25bp.bed -fo ${sample_size}_uncut_ext25bp.fa

	rm -f ${base}_random.bed

	#Handle on cut sites file
	bedtools slop -b 25 -i ${i} -g ${GENOME_SIZE_chrM} | bedtools getfasta -fi ${FASTA_FILE} -bed - \
	| awk -v RS=">" '!/N/{printf $0RT}' | grep '>' | sed 's/[:-]/\t/g' | sed 's/>//' \
	| awk -v FS="\t" -v OFS="\t" '{if(($3-$2) == 51) {print $0}}' | awk '!a[$1$2$3]++' \
	| shuf | head -${sample_size} | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,"1","."}' > ${sample_size}_all_ext25bp.bed
	bedtools getfasta -fi ${FASTA_FILE} -bed ${sample_size}_all_ext25bp.bed -fo ${sample_size}_all_ext25bp.fa

	#Merge cut and uncut files
	cat ${sample_size}_all_ext25bp.bed ${sample_size}_uncut_ext25bp.bed > ${base}_trainning_ext25bp.bed
	cat ${sample_size}_all_ext25bp.fa ${sample_size}_uncut_ext25bp.fa > ${base}_trainning_ext25bp.fa
	nohup python ~/Packages/BiasAway/BiasAway.py m -f ${base}_trainning_ext25bp.fa | sed "s/ >/>/" > ${base}_trainning_ext25bp_shuffled.fa &
done

ln -sf ~/scripts/generate_shapes.R .

for i in *_trainning_ext25bp.fa
do
	nohup Rscript generate_shapes.R -f ${i} &
done