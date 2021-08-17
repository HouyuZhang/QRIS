# **QRIS: a machine learning framework to investigate the determinants of retrovirus integration specificity**

**Author**: Houyu Zhang

**Email**: Hughiez047@gmail.com

#Copyright (c) 2021 YenLab@SKLEH. All rights reserved.

## Introduction

The activity of retrovirus is linked with various diseases. Retrovirus can integrate a copy of its genomic DNA into the host genome and store it as the provirus. It threatens the host cell by interrupting the host genome architecture and transcribing its provirus for detrimental expansion. Meanwhile, depending on where the retrovirus integrates, the host has corresponding defense mechanisms to repress the provirus transcription and hence eliminate the harm. So, the integration site selection is a vital process for retrovirus's fate. 

Many efforts were made to understand its integration specificity and found diverse DNA motif preferences across retroviruses genera. However, the effect of genome-wide DNA structural properties, like DNA shapes, on retrovirus integration was less clear. We systematically investigated this issue on six types of retroviruses that are representative of four genera. Here we devised QRIS (Quantify the Retrovirus Integration Specificity), a machine learning framework to assess the DNA shape effect on large-scale retroviruses integration sites. We found that the DNA shape can independently or cooperatively work with the DNA motif to regulate retrovirus integration. Based on this, we classified these retroviruses into three categories: StrongFavor retrovirus, WeakFavor retrovirus, and Strongshape retrovirus. Interestingly, the Strongshape retrovirus can gain specificity through DNA shapes even without a particular DNA motif. 

Our qualitative and quantitative evaluation of DNA shape and DNA motif revealed their diverse roles in regulating the retrovirus integration specificity. Our findings may help more precisely control the lentivirus vector for gene therapy and disturb the retrovirus integration during the pathogenic process in the future. 

![Graphical abstract](E:\OneDrive\1_Yenlab\GitHub\QRIS\Graphicalabstract.png)

## Scripts organization

All raw, intermediate, final data and all scripts here can fully reproduce the paper.

- Scripts head with series number is for specific analysis follows the order in the retrovirus paper.

  Specifically, the `bash script` is for processing data and corresponding `R scripts` is used for downstream analysis, statistics and plotting.
  

- QRIS_Rawdata: Store raw data, shuffled data, and R object of DNA shape values. These files if for QRIS_Step1.
- QRIS_Results: store machine learning results and can be visualized using QRIS_Step2.

## Citation

QRIS: a machine learning framework to investigate the determinants of retrovirus integration specificity. Sci Sin Vitae  