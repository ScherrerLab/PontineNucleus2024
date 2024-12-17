# PontineNucleus2024

This is a repository for code used in the publication Chen et al. 2024 entitled "Neural circuit basis of placebo pain relief". 

# single-cell RNA-seq
  Code for scRNA-seq:
  * `DSSeurat_functions.R` - This file contains R functions used in the primary scRNA-seq analysis pipeline.
  * `PN_seq.R` - This file contains the complete R code to cluster cells and generate plots used for Chen et al. 2024 from raw cell x gene matrices.

Raw fastq files and cells x genes matrices can be found at the following Gene Expression Omnibus (GEO) submission: GSE267264 for SMARTseq data and GSE267265 for 10X Genomics data.

# Ca2+ imaging 
Ca2+ activity of neurons extracted
* `rACC_Pn data` - This file contains the Ca2+ imaging data from rACC-Pn neurons procssed with EXTRACT.
* `PC_data` - This file contains the Ca2+ imaging data from PC neurons processed with EXTRACT.

# Custom R code
R code for Ca2+ iamging, behavior, and ephys data analysis
* `Placebo_state_final.R` - This file contains code for analysing rACC-Pn neuron Ca2+ imaging data.
* `PN_CB_PC_cluster.R` - This file contains code for analysing PC data.
* `PN_behaviro_analysis.R` - This file contains code for behavioral data analysis.
* `PN_ephys_04122021.R` - This file contains code for ephys data analysis.


#
Direct any questions regarding scRNA-seq analysis on this repository or email jesse_niehaus@med.unc.edu, any other questions to chong_chen@med.unc.edu
