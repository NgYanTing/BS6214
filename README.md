# BS6214 Code records for learning diaries

**(1) RNA-sequencing of mouse brains**

This was an attempt at analysing RNA-sequencing data from this [paper](https://doi.org/10.1016/j.celrep.2018.10.070) which reported the effects of various diet conditions on brain health. The task was to replicate the figures in the paper using a different set of analysis methods to see if the results would be different. Codes include processing of the data from 30 samples aligned to the mouse genome build mm39 with Kallisto-0.46.1. Analyses include differential gene expression analysis and GSEA using clusterProfiler.

It was found that the results related to RNA-sequencing reported in the paper were incredibly robust!


**(2) Population genetics of _Vitis vinifera_**

This was an exercise using a small subset of RNA-sequencing data of wine grapes (<i>Vitis vinifera</i>) from this [paper](https://www.nature.com/articles/s41467-021-27487-y). 8 wild grape varieties were compared to 11 cultivars to look for signs of selection on chromosome 17 as reported in the paper. Codes here are the R codes processing the outputs from GATK-4.1.9.0 (joint genotyping), PLINK-1.9 (PCA and IBS distances) and vcftools-0.1.16 (nucleotide diversities and fixation indices).

Analyses identified a similar region on chromosome 17 as reported in the paper that had likely undergone a selective sweep.
