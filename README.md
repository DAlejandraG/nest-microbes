# nest-microbes
This repository contains the code and data files needed to reproduce the analysis from Social organization and physical environment shape the microbiome of harvester ants" by Denisse Alejandra Gamboa, Peter Flynn, and Noa Pinter-Wollman.

Once the code and data are saved to your computer, you will need to change the file paths in the script to match your local directory. The lines that need editing should be obvious - they are the lines that call read.csv on a filepath through "C:/Users/Dalegamboa/Desktop/Projects/...". Once those read.csv function calls are changed to reflect the location of the data files on your computer, the R script in this repo should run without further modification. You'll see that the script begins with loading in data files and cleaning data objects, then proceeds to the analyses performed in the study, and concludes with the outputs found in the main text and supplementary material.

CleanUp_239.R is all R code used for data cleaning
Nest_microbes_code.R is all R code used in this project, including statistics and figures.
QIIME2_analysis.R is the QIIME2 code and functions used for this project. Please note that this file needs to be run in the command line for qiime2, not in R.

The data files:
sample_metadata_diversity_metrics_264.csv - contains information about the sample type, the nest ID, depth (cm), and values for diversity measures. There are the samples before data cleaning.

Data_Digging - sorted and marked for processing prioritization - with chamber type info.csv - contains information about the chamber types.

updated_clean_sample_metadata_diversity_metrics_239_1.csv - Post-cleaning results. Metadata for final 239 samples that will be used in the analysis. 

