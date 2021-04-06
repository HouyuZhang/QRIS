#!/bin/bash
##########################################################################
#This script is for retroviruse local pattern analysis
#Issue report on Hughiez047@gmail.com
#Copyright (c) 2020 __YenLab@SKLEH__. All rights reserved.
##########################################################################
#Usage: bash analyze_retroviruse.sh mm10 

source /public/home/zhy/scripts/0.utilities.sh
ln -sf /public/home/zhy/scripts/generate_shapes.R .

for i in arch/*.bed
do
	#Make sure the scafford and patch are removed from bed file
	#grep -vE ${CONSENSUS_REMOVE} ${i} > tmp && mv tmp ${i}

	base=`basename ${i} .bed`
	#Get species information for current sample
	id=`echo ${base} | sed "s/_.*//"`
	Configuration_info ${id} 
	echo "--->>> Start Process ${i} for ${GENOME_NAME} ..."

	#Make 51bp file
	ran1=25
	shuf ${i} | head -200000 | sort -k1,1 -k2,2n | bedtools slop -b ${ran1} -i - -g ${GENOME_SIZE_chrM} | \
	awk -v FS="\t" -v OFS="\t" '{if($3-$2 == 51) print}' > ${base}_ext${ran1}bp.bed
	bedtools getfasta -fi ${FASTA_FILE} -bed ${base}_ext${ran1}bp.bed > ${base}_ext${ran1}bp.fa

	nohup weblogo --composition ${COMPOSITION} -F pdf -S 0.2 --errorbars NO --color '#109648' A 'A' --color '#255C99' C 'C' --color '#D62839' T 'T' --color '#F7B32B' G 'G' \
    < ${base}_ext${ran1}bp.fa > ${base}_ext${ran1}bp_weblogo.pdf --stacks-per-line 51 &
    nohup meme -mod zoops -bfile ${markov_bg1} -p 20 ${base}_ext${ran1}bp.fa -oc ${base}_ext${ran1}bp_meme -revcomp -w 51 -searchsize 5000000 -nmotifs 3 -seed 0616 -dna &

	#Make 201bp file
	ran2=100
	shuf ${i} | head -200000 | sort -k1,1 -k2,2n | bedtools slop -b ${ran2} -i - -g ${GENOME_SIZE_chrM} |\
	awk -v FS="\t" -v OFS="\t" '{if($3-$2 == 201) print}' > ${base}_ext${ran2}bp.bed
	bedtools getfasta -fi ${FASTA_FILE} -bed ${base}_ext${ran2}bp.bed > ${base}_ext${ran2}bp.fa
	nohup Rscript generate_shapes.R -f ${base}_ext${ran2}bp.fa &
	#Rscript ~/Zhang_Scripts/Zhang/R/Fasta_GC_content.R -f ${base}_ext${ran2}bp.fa -s 1

	if [ ! -f ${GENOME_NAME}.random_ext25bp.fa ]
	then
		bedtools random -l 51 -n 200000 -g ${GENOME_SIZE_chrM} > ${GENOME_NAME}.random_ext25bp.bed
		bedtools getfasta -fi ${FASTA_FILE} -bed ${GENOME_NAME}.random_ext25bp.bed > ${GENOME_NAME}.random_ext25bp.fa

		bedtools random -l 1001 -n 200000 -g ${GENOME_SIZE_chrM} > ${GENOME_NAME}.random_ext${ran2}bp.bed
		bedtools getfasta -fi ${FASTA_FILE} -bed ${GENOME_NAME}.random_ext${ran2}bp.bed > ${GENOME_NAME}.random_ext${ran2}bp.fa
		nohup Rscript generate_shapes.R -f ${GENOME_NAME}.random_ext${ran2}bp.fa &
		#Rscript ~/Zhang_Scripts/Zhang/R/Fasta_GC_content.R -f ${GENOME_NAME}.random_ext${ran2}bp.fa -s 1
	fi
done

zip -r all_meme_result.zip *meme*
rm -f *fa.*

for i in *ext25bp.bed
do
	#base=`basename ${i} .bed`
	#GENOME_SIZE_chrM="/public/home/zhy/Genomes_info/homo_sapiens/hg38.chrom.sizes"
	#FASTA_FILE="/public/home/zhy/Genomes_info/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly_chrM.fa"
	#COMPOSITION=38.84%
#
	#GENOME_SIZE_chrM="/public/home/zhy/Genomes_info/Mus_musculus/mm10.chrom.sizes"
	#FASTA_FILE="/public/home/zhy/Genomes_info/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa"
	#COMPOSITION=40.48%

	GENOME_SIZE_chrM="/public/home/zhy/Genomes_info/danio_rerio/danRer11.chrom.sizes"
	FASTA_FILE="/public/home/zhy/Genomes_info/danio_rerio/Danio_rerio.GRCz11.dna.primary_assembly_chrM.fa"
	COMPOSITION=36.52%

	awk -v FS="\t" -v OFS="\t" '{if($3-$2 == 51) print}' ${i} | bedtools getfasta -fi ${FASTA_FILE} -bed - > ${base}.fa
    #meme -mod zoops -bfile ~/Genomes_info/homo_sapiens/hg38_markov_bg_1 -p 20 ${i} -oc ${base}_meme -revcomp -w 51 -searchsize 5000000 -nmotifs 3 -seed 0616 -dna 
    weblogo --composition ${COMPOSITION} -F pdf -S 0.2 --errorbars NO --color '#109648' A 'A' --color '#255C99' C 'C' --color '#D62839' T 'T' --color '#F7B32B' G 'G' \
    -f ${base}.fa -o ${base}_weblogo.pdf --stacks-per-line 51
done










