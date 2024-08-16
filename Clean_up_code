### Data clean up for "Social organization and physical environment shape the microbiome of harvester ants"
## authors: Alejandra Gamboa, Peter Flynn, Noa Pinter-Wollman
## Email Alejandra Gamboa at dgamboa24@g.ucla.edu with questions on this script and data files

#################################### START UP ###################################################### 
# clear workspace, load packages
rm(list = ls())
# set up
library(stringr)
library(dplyr)
library(tidyr)
# load data
metadata = read.csv("sample_metadata_diversity_metrics_264.csv")
chamberIDdata = read.csv("Data_Digging - sorted and marked for processing prioritization - with chamber type info.csv", header = T)

##### clean up data entry errors per lab notes from lab notebook, extraction sheets, and field data notes 
# rename C-G2-1 to C-G3 for col 'id'
fix_names_CG3 = metadata
fix_names_CG3$id[fix_names_CG3$id == 'C-G2-1'] <- 'C-G3'
# Replace for col:'chamber-Iden' 'sample.ID.on.tube' by checking condition on fixed col 'id'
fix_names_CG3$Chamber_Iden[fix_names_CG3$id == 'C-G3'] <- "G3"
fix_names_CG3$sample.ID.on.tube[fix_names_CG3$id == 'C-G3'] <- "C G3"

# rename E-I3 to E-J1 (samples swapped)
fix_names_EJ1 = fix_names_CG3
fix_names_EJ1$id[fix_names_EJ1$id == 'E-I3'] <- 'E-J1'
# Replaced col by checking condition on fixed col 'id'
fix_names_EJ1$Chamber.ID[fix_names_EJ1$id == 'E-J1'] <- "J" # correct Chamber.ID
fix_names_EJ1$Chamber_Iden[fix_names_EJ1$id == 'E-J1'] <- "J1" # correct chamber_Iden 
fix_names_EJ1$sample.ID.on.tube[fix_names_EJ1$id == 'E-J1'] <- "E J1" # correct sample.ID.on.tube
fix_names_EJ1$Samples[fix_names_EJ1$id == 'E-J1'] <- "ants" # correct Samples
fix_names_EJ1$Sample_Type[fix_names_EJ1$id == 'E-J1'] <- "ants" # correct Samples_Type
fix_names_EJ1$Sample_Type2[fix_names_EJ1$id == 'E-J1'] <- "ants"# correct Sample_Type2
fix_names_EJ1$Nest_Sample[fix_names_EJ1$id == 'E-J1'] <- "Eggplant_ants"# correct Nest_Sample
fix_names_EJ1$Depth_categorical[fix_names_EJ1$id == 'E-J1'] <- "medium" # correct Depth_categorical
fix_names_EJ1$depth_nest[fix_names_EJ1$id == 'E-J1'] <- "Eggplant_27" # correct depth_nest
fix_names_EJ1$depth_categorical[fix_names_EJ1$id == 'E-J1'] <- 27 # correct depth_categorical
fix_names_EJ1$depth_continous[fix_names_EJ1$id == 'E-J1'] <- 0.54 # correct depth_continous
fix_names_EJ1$sample_depth[fix_names_EJ1$id == 'E-J1'] <- "ants27" # correct sample_depth
fix_names_EJ1$chamberid_depth[fix_names_EJ1$id == 'E-J1'] <- "J_27" # correct chamberid_depth
fix_names_EJ1$samples_good[fix_names_EJ1$id == 'E-J1'] <- "ants"# correct samples_good
fix_names_EJ1$nest[fix_names_EJ1$id == 'E-J1'] <- "27" # correct nest

# rename E-I3-1-to E-I3 (samples swapped; duplicate names when sent for sequencing)
fix_names_EI3 = fix_names_EJ1
fix_names_EI3$id[fix_names_EI3$id == 'E-I3-1'] <- 'E-I3'
# Replaced col by checking condition on fixed col 'id'
fix_names_EI3$Samples[fix_names_EI3$id == 'E-I3'] <- "soil" # correct Samples
fix_names_EI3$Sample_Type[fix_names_EI3$id == 'E-I3'] <- "soil" # correct Samples_Type
fix_names_EI3$Sample_Type2[fix_names_EI3$id == 'E-I3'] <- "soil"# correct Sample_Type2
fix_names_EI3$Nest_Sample[fix_names_EI3$id == 'E-I3'] <- "Eggplant_soil"# correct Nest_Sample 
fix_names_EI3$sample_depth[fix_names_EI3$id == 'E-I3'] <- "soil12"# correct sample_depth 
fix_names_EI3$samples_good[fix_names_EI3$id == 'E-I3'] <- "soil"# correct samples_good

# fix depth for E-I5 
fix_depth_EI5 = fix_names_EI3
fix_depth_EI5$Depth_categorical[fix_depth_EI5$id == 'E-I5'] <- "high" # correct Depth_categorical
fix_depth_EI5$depth_nest[fix_depth_EI5$id == 'E-I5'] <- "Eggplant_12" # correct depth_nest
fix_depth_EI5$depth_categorical[fix_depth_EI5$id == 'E-I5'] <- 12 # correct depth_categorical
fix_depth_EI5$depth_continous[fix_depth_EI5$id == 'E-I5'] <- 0.24 # correct depth_continous
fix_depth_EI5$sample_depth[fix_depth_EI5$id == 'E-I5'] <- "control far12" # correct sample_depth
fix_depth_EI5$chamberid_depth[fix_depth_EI5$id == 'E-I5'] <- "I_12" # correct chamberid_depth
fix_depth_EI5$nest[fix_depth_EI5$id == 'E-I5'] <- "12" # correct nest

# rename DC-F2-1 to DC-F3
fix_names_DCF3 = fix_depth_EI5
fix_names_DCF3$id[fix_names_DCF3$id == 'DC-F2-1'] <- 'DC-F3'
# Replace for col:chamber-Iden' 'sample.ID.on.tube' by checking condition on fixed col 'id'
fix_names_DCF3$Chamber_Iden[fix_names_DCF3$id == 'DC-F3'] <- "F3"
fix_names_DCF3$sample.ID.on.tube[fix_names_DCF3$id == 'DC-F3'] <- "DC F3"

# rename E-M2-1 to E-M2 
fix_names_EM2 = fix_names_DCF3
fix_names_EM2$id[fix_names_EM2$id == 'E-M2-1'] <- 'E-M2'
# Replace for col:'chamber-Iden' 'sample.ID.on.tube' by checking condition on fixed col 'id'
fix_names_EM2$Chamber_Iden[fix_names_EM2$id == 'E-M2'] <- "M2"
fix_names_EM2$sample.ID.on.tube[fix_names_EM2$id == 'E-M2'] <- "E M2"
###################### E-D1 IS CORRECT NO CHANGE IT IS BROOD
###################### E-D1-1 IS SOIL E-D2
# rename E-D1-1 to E-D2  (E-D2 SOIL sample; there were 2 extractions for E-D2, soil and ants)
fix_names_ED2 = fix_names_EM2
fix_names_ED2$id[fix_names_ED2$id == 'E-D1-1'] <- 'E-D2'
# Replace col by checking condition on fixed col 'id'
fix_names_ED2$Chamber_Iden[fix_names_ED2$id == 'E-D2'] <- "D2" # correct Chamber_Iden
fix_names_ED2$sample.ID.on.tube[fix_names_ED2$id == 'E-D2'] <- "E D2" # correct sample.ID.on.tube
fix_names_ED2$Samples[fix_names_ED2$id == 'E-D2'] <- "soil" # correct Samples
fix_names_ED2$Sample_Type[fix_names_ED2$id == 'E-D2'] <- "soil" # correct Sample_Type
fix_names_ED2$Nest_Sample[fix_names_ED2$id == 'E-D2'] <- "Eggplant_soil" # correct Nest_Sample
fix_names_ED2$sample_depth[fix_names_ED2$id == 'E-D2'] <- "soil28" # correct sample_depth

# change name of 'depth_categorical' to 'depth' because column was mislabeled
corr_depth = fix_names_ED2
colnames(corr_depth)[which(names(corr_depth) == "Depth")] <- "Depth_categorical"
colnames(corr_depth)[which(names(corr_depth) == "depth_categorical")] <- "Depth"

# correct row, replace 'leaf' with 'seeds'
cor_seeds = corr_depth
cor_seeds$Samples = str_replace_all(cor_seeds$Samples, "leaf", "seeds")

# correct row, replace 'beetle/seed/food with 'seeds
cor_seeds2 = cor_seeds
cor_seeds2$Samples = str_replace_all(cor_seeds2$Samples, "beetle/seeds/food", "seeds")

# remove extra empty column
rm_col = cor_seeds2[,-27] # empty col

# clean up by 'Samples', remove 'entrance soil', 'soil control off mound", "soil control off mount", "soil control near mound" 
clean_meta = rm_col[!(rm_col$Samples=="entrance soil" | rm_col$Samples=="soil control off mound" | 
                        rm_col$Samples=="soil control near mound" | rm_col$Samples=="soil contol off mount"), ]

# remove samples with 'entrance chamber' type: Carrot J1, J2, J3, J4 (refer to spreadsheet with chamber type info) 
clean_meta2 = clean_meta[!((clean_meta$id) %in% c("C-J1","C-J2","C-J3","C-J4")),]

# remove samples with 'outside chamber' type (refer to spreadsheet with chamber type info) :
# Apple J1, J2, J3, J4, M1, M2, M3, M4, M5
# DC F1, F2, F3, F4, F5
# Eggplant M1, M2, M3, M4
clean_meta3 = clean_meta2[!((clean_meta2$id) %in% c("A-J1","A-J2","A-J3","A-J4","A-M1","A-M2","A-M3","A-M4","A-M5","DC-F1","DC-F2","DC-F3","DC-F4","DC-F5","E-M1","E-M2","E-M3","E-M4")),]

# remove 3 rows for Banana A1, A2, A3 because they are not real chamber samples
# remove sample A-H2 because this sample contained a mix of soil and seeds
clean_meta4 = clean_meta3[!((clean_meta3$id) %in% c("B-A1","B-A2","B-A3","A-H2")),]

#### Prepare to merge data frames
# 1. Clean up chamberIDdata before merging with clean meta data
# remove extra spaces from col 'Chamber.Type.simple 
chamberIDdata$Chamber.Type.simple = str_trim(chamberIDdata$Chamber.Type.simple)
# replace blank strings with NA in 'chamber.type.simple' col
#Replace only character columns
clean_chamberID_data = chamberIDdata
clean_chamberID_data$Chamber.Type.simple[clean_chamberID_data$Chamber.Type.simple == ''] <- NA 

# 2. Combine "Nest" and "ChamberIden" columns in "clean_meta4" df to get a "Sample.ID" column that matches "clean_chamberID_data" 
clean_metadata = clean_meta4
clean_metadata$Sample.ID = paste(clean_metadata$Nest, clean_metadata$Chamber_Iden)

### Merge "Chamber.Type.Simple" column to "clean_metadata" by "Sample.ID" column
mergedMetadata= merge(clean_metadata, clean_chamberID_data[ ,c("Sample.ID", "Chamber.Type.simple")], by= "Sample.ID", all.x = TRUE)

## fix 'DC-Q4'chamber type simple is 'brood/reproductive' NOT 'NA
merged_metadata_clean = mergedMetadata
merged_metadata_clean$Chamber.Type.simple[merged_metadata_clean$id == 'DC-Q4'] <- "brood/reproductive"

# export clean data as csv 
write.csv(merged_metadata_clean, file = "updated_clean_sample_metadata_diversity_metrics_239.csv")
