# nest-microbes
This repository contains the code and data files needed to reproduce the analysis from "Social organization and physical environment shape the microbiome of harvester ants" by Denisse Alejandra Gamboa, Peter Flynn, and Noa Pinter-Wollman.

Once the code and data are saved to your computer, you will need to change the file paths in the script to match your local directory. Once the read.csv function calls are changed to reflect the location of the data files on your computer, the R script in this repo should run without further modification. You'll see that the script begins with loading in data files and cleaning data objects, then proceeds to the analyses performed in the study, and concludes with the outputs found in the main text and supplementary material.

**Code:**
1. QIIME2_analysis.R is the QIIME2 code and functions used for this project. Please note that this file needs to be run in the command line for QIIME2, not in R.
This code takes sequencing data from the NCBI reposetory found at NCBI BioProject record SUB14655354, https://www.ncbi.nlm.nih.gov/sra/PRJNA1147938 and outputs summary statistics of microbiome diversity, including the various alpha and beta diversity measures.

2. Nest_microbes_code.R is R code used to analyse data and produce figures for the manuscript. It requires the file "updated_clean_sample_metadata_diversity_metrics_239_1.csv"

**Data files (in folder 'Data':**
1. table.qza - Amplicon Aequence Variant (ASV) table, which records the number of times each exact amplicon sequence variant was observed in each sample. It is produced by the "QIIME2_analysis.R" code and is required for the "Nest_microbes_code.R".
   
2. taxonomy.qza - Contains taxonomic assignments data. It is produced by the "QIIME2_analysis.R" code and is required for the "Nest_microbes_code.R".
   
3. rooted-tree.qza - File that contains a phylogenetic tree that has been rooted and is important for phylogenetic diversity computation, such as Faith's Phylogenetic Diversity and UniFraq. It is produced by the "QIIME2_analysis.R" code and is required for the "Nest_microbes_code.R".

4. updated_clean_sample_metadata_diversity_metrics_239_1.csv - Metadata for the 239 samples included in the data that is analyzed by "Nest_microbes_code.R".

5. data_descriptions - a description of all the columns in "updated_clean_sample_metadata_diversity_metrics_239_1.csv"

