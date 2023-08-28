# Workflow 
This page decribes scripts and data used for the analysis of PTG (parathyroid gland) chromatin and expression data. 

Due to the data size, datasets required for `1_expression_gcm2_analysis.r` and `4_ptg_genes_chromatin_expression.r` can be downloaded [here](https://www.dropbox.com/scl/fo/rtwdt0nd6e6p1jg666gph/h?rlkey=cesecj7intlh9zc1y9czlnumz&dl=0).

## PTG expression compared to other tissues
`1_expression_gcm2_analysis.r` takes input of gene expression levels from 8 PTG individuals (RNA-seq) and compared them with those from other tissues in GTEX. `f1.data.RData` contains PTG gene expression data as well GTEX data. The clustering and visualization was done based on PTG-spepcifically expressed genes. 

## PTG chromatin interactions compared to those in other tissues
`2_chromatin_interaction_analysis.r` takes input of chromatin interactions of Hi-C data from 2 PTG individuals and compared them with those from other tissues (pcHi-C) in [Jung et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6778519/). Significant interactions were detected using [Juicer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5846465/) HiCCUPS with FDR = 0.1 and 5 kb resolution for the merged data from biological replicates.
