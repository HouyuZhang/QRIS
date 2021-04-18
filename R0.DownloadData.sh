#!/bin/bash

nohup curl -L https://rid.ncifcrf.gov/tmp/1101457718.dwnld.txt -o HIV.txt &
nohup curl -L https://rid.ncifcrf.gov/tmp/2097771456.dwnld.txt -o HTLV1.txt &
nohup curl -L https://rid.ncifcrf.gov/tmp/1529399652.dwnld.txt -o SIV.txt &

cut -f 3-4 HTLV1.txt | tail -n +2 | awk -v FS="\t" -v OFS="\t" '{{print $1,$2-1,$2}}' \
| awk '!a[$1$2$3]++' | sort -k1,1 -k2,2n > Human_HTLV-1_Tcells_2017Artesi.bed
liftOver Human_HTLV-1_Tcells_2017Artesi.bed ~/Genomes_info/homo_sapiens/hg19ToHg38.over.chain.gz tmp unlifted.bed
mv tmp Human_HTLV-1_Tcells_2017Artesi.bed

cut -f 3-4 SIV.txt | tail -n +2 | awk -v FS="\t" -v OFS="\t" '{{print $1,$2-1,$2}}' \
| awk '!a[$1$2$3]++' | sort -k1,1 -k2,2n > Human_SIV_2019Ferris.bed
liftOver Human_SIV_2019Ferris.bed ~/Genomes_info/homo_sapiens/hg19ToHg38.over.chain.gz tmp unlifted.bed
mv tmp Human_SIV_2019Ferris.bed

awk -v FS="\t" -v OFS="\t" '{print $0 >> $17".bed"}' HIV.txt
for i in *bed
do
	cut -f 3-4 ${i} | awk -v FS="\t" -v OFS="\t" '{{print $1,$2-1,$2}}' \
	| awk '!a[$1$2$3]++' | sort -k1,1 -k2,2n > tmp
	liftOver tmp ~/Genomes_info/homo_sapiens/hg19ToHg38.over.chain.gz ${i}_ unlifted.bed
done

#=========================================================================================
# Handle dowloaded file
#=========================================================================================
#1. HIV file 
#(Wang et al., 2007)
liftOver hiv_Avr.wig.bed ~/Genomes_info/homo_sapiens/hg19ToHg38.over.chain.gz HIV-1_Jurkat_Avr_2007Wang.bed unlifted.bed
liftOver hiv_Mse.wig.bed ~/Genomes_info/homo_sapiens/hg19ToHg38.over.chain.gz HIV-1_Jurkat_Mse_2007Wang.bed unlifted.bed

grep -vE "GL|JH|KI|GL|chrPt|KN|KZ|Scaffold|chr2110|chrrDNA|ctg" HIV-1_Jurkat_Avr_2007Wang.bed >tmp && mv tmp HIV-1_Jurkat_Avr_2007Wang.bed
grep -vE "GL|JH|KI|GL|chrPt|KN|KZ|Scaffold|chr2110|chrrDNA|ctg"  HIV-1_Jurkat_Mse_2007Wang.bed  >tmp && mv tmp  HIV-1_Jurkat_Mse_2007Wang.bed

#(Sherrill-Mix et al., 2013) Active, Bcl-2 transduced, Central Memory,Jurkat, Resting types of T cell
cut -d "," -f 2-5 12977_2013_3591_moesm2_esm | awk -v FS="," -v OFS="\t" '{print $1,$2-1,$2,$4,".",$3}' |\
 sed "s/\"//g" |sed -e '1d' > HIV-1_Tcells_2013Sherrill.bed

#2. MLV file
#(LaFave et al., 2014)
cut -f 1-6 mlv_k562_LaFave_et_al.bed | sed -e '1d' > tmp 
liftOver tmp ~/Genomes_info/homo_sapiens/hg19ToHg38.over.chain.gz MLV_K562_2014LaFave.bed unlifted.bed
cut -f 1-6 mlv_hepg2_LaFave_et_al.bed | sed -e '1d' > tmp  
liftOver tmp ~/Genomes_info/homo_sapiens/hg19ToHg38.over.chain.gz MLV_Hep2_2014LaFave.bed unlifted.bed

grep -vE "GL|JH|KI|GL|chrPt|KN|KZ|Scaffold|chr2110|chrrDNA|ctg" MLV_K562_2014LaFave.bed | awk -v FS="\t" -v OFS="\t" '{print $1,$2+2,$3-1,$4,$5,$6}'|\
 awk '!a[$1$2$3$6]++' > tmp && mv tmp MLV_K562_2014LaFave.bed
grep -vE "GL|JH|KI|GL|chrPt|KN|KZ|Scaffold|chr2110|chrrDNA|ctg" MLV_Hep2.bed | awk -v FS="\t" -v OFS="\t" '{print $1,$2+2,$3-1,$4,$5,$6}'|\
 awk '!a[$1$2$3$6]++' > tmp && mv tmp MLV_Hep2_2014LaFave.bed

#3. HTLV-1
#(Furuta et al., 2017) from a word file
awk -v FS="\t" -v OFS="\t" '{print $1,$2-1,$2,".",".",$3}' data.txt > HTLV-1_hematopoietic_2017Furuta.bed 

#4. MMTV, SB, PB
#de Jong et al., 2014
for i in MouseMammaryTumorVirus_NMuMG_2014Jong.bed SleepingBeauty_mESC_2014Jong.bed piggyBac_mESC_2014Jong.bed;\
do awk '!a[$1$2$3$6]++' ${i} > tmp && mv tmp ${i};done

#5. MMLV
#(Varshney et al. 2013) danRer7 -> danRer10 -> danRer11
liftOver 2013_07_01_zebrafish_insertion_bedfiles.bed ~/Genomes_info/danio_rerio/danRer7ToDanRer10.over.chain.gz tmp unlifted.bed
liftOver tmp ~/Genomes_info/danio_rerio/danRer10ToDanRer11.over.chain.gz MMLV_Blastula_2013Varshney.bed unlifted.bed
grep -vE "GL|JH|KI|GL|chrPt|KN|KZ|Scaffold|chr2110|chrrDNA|ctg" MMLV_Blastula_2013Varshney.bed >tmp && mv tmp MMLV_Blastula_2013Varshney.bed

# VISDB (JTLV-1, HIV) Many Wrong data!!!
for i in *.xlsx; do libreoffice --headless --convert-to csv ${i} ; done

awk -v FS="," -v OFS="\t" '{print $6,$9,$10,$11,$13,$14 >> "VIS-HIV_"$14".bed"}' VIS-HIV.csv
awk -v FS="," -v OFS="\t" '{print $6,$9,$10,$11,$13,$14 >> "VIS-HTLV-1_"$14".bed"}' VIS-HTLV-1.csv
