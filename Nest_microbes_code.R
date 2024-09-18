### Analysis for "Social organization and physical environment shape the microbiome of harvester ants"
## authors: Alejandra Gamboa, Peter Flynn, Noa Pinter-Wollman
## Email Alejandra Gamboa at dgamboa24@g.ucla.edu with questions on this script and data files

#################################### START UP ###################################################### 
# clear workspace, load packages
rm(list = ls())
## stats
library('performance')
library('ggplot2')
library('car')
library("emmeans")
library("vegan")
library("pairwiseAdonis")
library("multcomp")
## FIGURES in R 
library("qiime2R")
library("phyloseq")
library("microbiome")
library("microbiomeutilities")
library("ggtext")
library("ggraph")
library("DT")
library("corncob")
library("microViz")
library("devtools")
library("ggpubr")
library("ggplot2")
library("gridExtra")
library("patchwork")

############################# Load data files, set working directory, create/edit objects needed:
# load metadata file 
metadata <- read.csv("updated_clean_sample_metadata_diversity_metrics_239_1.csv")
# convert Fig4_sample column in the metadata dataframe into a factor with the specified levels
metadata$Fig4_Sample = factor(metadata$Fig4_Sample, levels=c("ants", "reproductives", "brood", "seeds", "soil"))

# Import data from QIIME 2 artifacts into a phyloseq object for microbial community/diversity analysis
phy <- qza_to_phyloseq("table.qza", "rooted-tree.qza", "taxonomy.qza","updated_clean_sample_metadata_diversity_metrics_239_1.tsv")

# create a theme for plotting that looks normal:
theme_mine <- function(base_size = 18, base_family = "Helvetica") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 25),
      strip.text.y = element_text(size = 25),
      axis.text.x = element_text(size=25),
      axis.text.y = element_text(size=25,hjust=1),
      axis.ticks =  element_line(colour = "black"),
      axis.title.x= element_text(size=25),
      axis.title.y= element_text(size=25,angle=90),
      panel.background = element_blank(),
      panel.border =element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.margin = unit(1.0, "lines"),
      plot.background = element_blank(),
      plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1)
    )
}

############### Analysis and figures for ALL SAMPLE TYPES #####################################
################################################################################################

##### FIGURE 3 - RELATIVE ABUNDANCE #######
pseq <-  phy %>%  # fix taxonomy issue and validate the fixed phyloseq object
  tax_fix() %>%
  phyloseq_validate() 
pseq <- pseq %>% tax_fix(unknowns = c("uncultured")) 

# Create ordered factor for `samples_reorder` in the sample data of `pseq` and `phy`
sample_levels <- c("ants", "reproductives", "brood", "seeds", "soil", "control near", "control far")
pseq@sam_data$samples_reorder <- factor(pseq@sam_data$Samples, levels = sample_levels)
phy@sam_data$samples_reorder <- factor(phy@sam_data$Samples, levels = sample_levels)

# make new labels and colors for barplot
sample_labels = c("ants" = "Ants","reproductives" = "Rep","brood" = "Brood","seeds" = "Seeds",
  "soil" = "Nest soil","control near" = "Control Near","control far" = "Control Far")

sample_colors = c("ants" = "darkviolet","reproductives" = "plum1","brood" = "mediumorchid1",
  "seeds" = "olivedrab3", "soil" = "goldenrod3","control near" = "bisque3","control far" = "bisque1")

# define a custom labeller function to rename and color background labels
custom_labeller <- function(variable, value) {
  # Define a custom labeller function
  value <- sample_labels[value]
  return(value)
}

phy_rel_total <- pseq %>% # create comp bar plot
  comp_barplot(
    tax_level = "Phylum",
    label = "SAMPLE", # name an alternative variable to label axis
    n_taxa = 32, # give more taxa unique colours
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other phylum", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  )  +  theme(axis.text.x=element_text(angle=90)) + theme(legend.position = "right") + facet_grid(~samples_reorder, scales = "free", space = "free", labeller = labeller(samples_reorder = custom_labeller))

phy_rel_total= phy_rel_total + labs(
  y = "Relative Abundance", x = "Samples", 
    ) + theme(
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 15))

##### run models - LM: diversity ~ sample type (use Fig4_Sample column) ####
############################################################################
# shannon
lmShannon = lm(shannon ~ Fig4_Sample, data = metadata)
check_model(lmShannon)
summary(lmShannon)
Anova(lmShannon)

# faith_pd
lmFaith = lm(faith_pd ~ Fig4_Sample, data = metadata)
check_model(lmFaith)
summary(lmFaith)
Anova(lmFaith)

# pielou_evenness
lmPE = lm(pielou_evenness ~ Fig4_Sample, data = metadata)
check_model(lmPE)
summary(lmPE)
Anova(lmPE)

# observed
lmObserved = lm(observed_otus ~ Fig4_Sample, data = metadata)
check_model(lmObserved)
summary(lmObserved)
Anova(lmObserved)
 
########## post hoc tukey test #############################
#################################################################
aovShannon = aov(shannon ~ Fig4_Sample, data = metadata)
check_model(aovShannon)
summary(aovShannon)
TukeyHSD(aovShannon)

aovFaith = aov(faith_pd ~ Fig4_Sample, data = metadata)
check_model(aovFaith)
summary(aovFaith)
TukeyHSD(aovFaith)

aovPE = aov(pielou_evenness ~ Fig4_Sample, data = metadata)
check_model(aovPE)
summary(aovPE)
TukeyHSD(aovPE)

aovObserved = aov(observed_otus ~ Fig4_Sample, data = metadata)
check_model(aovObserved)
summary(aovObserved)
TukeyHSD(aovObserved)

############ PERMANOVA for ord plot ######################
###########################################################
# Extract sample data and distance matrix
sample_data1 = as(sample_data(pseq), "data.frame")
bray_dist1 = distance(pseq, method = "bray")

# Perform PERMANOVA
permanova_result1 = adonis2(bray_dist1 ~ Fig4_Sample, data = sample_data1, permutations = 999)
print(permanova_result1)

# perform pairwise permanova tests
pairwise_result1 = pairwise.adonis(bray_dist1, sample_data1$Fig4_Sample)
print(pairwise_result1)

##########################################################################
############ FIGURE 4 - Diversity Metrics and Overall PCOA ################
#################### SHANNON PLOT #################
# Data frame for the post-hoc labels
max_value = max(metadata$shannon) # Determine the maximum value of the y-axis to place the labels uniformly
label_y_position = max_value + 1  
labels = data.frame(
  Fig4_Sample = c("ants", "reproductives", "brood", "seeds", "soil"),
  shannon = rep(label_y_position, 5),  # Use the same y position for all labels 
  label = c("A", "B", "B", "C", "D")
)
Shannon <- ggboxplot(metadata, x = "Fig4_Sample", y = "shannon", color = "Fig4_Sample") + scale_x_discrete(labels = c("Ants", "Repro", "Brood", "Seeds", "Soil")) + ggtitle("Shannon")
Shannon = Shannon + geom_boxplot(fill = c( "darkviolet", "plum1", "mediumorchid1", "olivedrab3", "sienna")) + theme_mine() + ggtitle("Shannon") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 35), legend.position = "none") + theme(axis.title.x = element_blank()) +theme(axis.text.x=element_text(size=30, angle = 0), axis.text.y = element_text(size = 25, angle = 90), axis.title.y = element_text(size = 30)) + 
  geom_text(data = labels, aes(x = Fig4_Sample, y = shannon, label = label), 
            vjust = -0.5, size = 15) + ylim(NA, label_y_position + 1) + # set y-axis limits 
            ylab("Shannon")  
  Shannon

#################### FAITH PLOT #################
# Data frame for the post-hoc labels
max_value = max(metadata$faith_pd) # Determine the maximum value of the y-axis to place the labels uniformly
label_y_position = max_value + 1  
labels = data.frame(
  Fig4_Sample = c("ants", "reproductives", "brood", "seeds", "soil"),
  faith_pd = rep(label_y_position, 5),  # Use the same y position for all labels 
  label = c("A", "B", "BC", "D", "C")
)
faith_pd <- ggboxplot(metadata, x = "Fig4_Sample", y = "faith_pd", color = "Fig4_Sample") + scale_x_discrete(labels = c("Ants", "Repro", "Brood", "Seeds", "Soil")) + ggtitle("Faith PD")
faith_pd = faith_pd + geom_boxplot(fill = c( "darkviolet", "plum1", "mediumorchid1", "olivedrab3", "sienna")) + theme_mine() + ggtitle("Faith PD") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 35), legend.position = "none") + theme(axis.title.x = element_blank()) +theme(axis.text.x=element_text(size=30, angle = 0),axis.text.y = element_text(size = 25, angle = 90), axis.title.y = element_text(size = 30)) + 
  geom_text(data = labels, aes(x = Fig4_Sample, y = faith_pd, label = label), 
            vjust = -0.5, size = 15) + ylim(NA, label_y_position + 2) + # set y-axis limits 
            ylab("Faith's PD")
faith_pd

#################### EVENNESS PLOT #################
# Data frame for the post-hoc labels
max_value = max(metadata$pielou_evenness) # Determine the maximum value of the y-axis to place the labels uniformly
label_y_position = max_value + 0.03  
labels = data.frame(
  Fig4_Sample = c("ants", "reproductives", "brood", "seeds", "soil"),
  pielou_evenness = rep(label_y_position, 5),  # Use the same y position for all labels 
  label = c("A", "BC", "B", "D", "C")
)
Even <- ggboxplot(metadata, x = "Fig4_Sample", y = "pielou_evenness", color = "Fig4_Sample") + scale_x_discrete(labels = c("Ants", "Repro", "Brood", "Seeds", "Soil")) + ggtitle("Evenness")
Even = Even + geom_boxplot(fill = c( "darkviolet", "plum1", "mediumorchid1", "olivedrab3", "sienna")) + theme_mine() + ggtitle("Evenness") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 35), legend.position = "none") + theme(axis.title.x = element_blank()) +theme(axis.text.x=element_text(size=30, angle = 0),axis.text.y = element_text(size = 25, angle = 90), axis.title.y = element_text(size = 30)) + 
  geom_text(data = labels, aes(x = Fig4_Sample, y = pielou_evenness, label = label), 
            vjust = -0.5, size = 15) + ylim(NA, label_y_position + 0.07) + # set y-axis limits
            ylab("Pielou's evenness")
Even

#################### OBSERVED #################
# Data frame for the post-hoc labels
max_value = max(metadata$observed_otus) # Determine the maximum value of the y-axis to place the labels uniformly
label_y_position = max_value + 25  
labels = data.frame(
  Fig4_Sample = c("ants", "reproductives", "brood", "seeds", "soil"),
  observed_otus = rep(label_y_position, 5),  # Use the same y position for all labels 
  label = c("A", "BC", "BD", "C", "D")
)
Observed <- ggboxplot(metadata, x = "Fig4_Sample", y = "observed_otus", color = "Fig4_Sample") + scale_x_discrete(labels = c("Ants", "Repro", "Brood", "Seeds", "Soil")) + ggtitle("Observed ASVs")
Observed = Observed + geom_boxplot(fill = c( "darkviolet", "plum1", "mediumorchid1", "olivedrab3", "sienna")) + theme_mine() + ggtitle("Observed ASVs") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 35), legend.position = "none") + theme(axis.title.x = element_blank()) +theme(axis.text.x=element_text(size=30, angle = 0),axis.text.y = element_text(size = 25, angle = 90), axis.title.y = element_text(size = 30)) + 
  geom_text(data = labels, aes(x = Fig4_Sample, y = observed_otus, label = label), 
            vjust = -0.5, size = 15) + ylim(NA, label_y_position + 2) + # set y-axis limits
          ylab("Observed ASVs")
Observed
# arrange all four plots into one image
quartz()
grid.arrange(Shannon, faith_pd, Even, Observed,  nrow= 2)
#pdf saved 16x11

################### FIGURE 4B total ord with repro - Bray-Curtis PCoA plot ########################
total_phy = pseq %>%
  tax_transform(rank = "unique", trans = "identity") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Fig4_Sample", fill = "Fig4_Sample",
    shape = "circle", alpha = 0.8,
    size = 2
  ) +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = Fig4_Sample)
  )
# add custom theme and color scales
total_phy <- total_phy + theme_mine() + scale_colour_manual(values = c( "darkviolet", "plum1", "mediumorchid1", "olivedrab3", "sienna")) + guides(colour = guide_legend(title = "Sample Type"), fill = guide_legend(title = "Sample Type")) + 
  theme(
  legend.position = c(0.7, 0.05), # adjust legend 
  legend.justification = c(0, 0), 
  legend.text = element_text(size = 30),
  axis.title.x = element_text(size = 30),
  axis.title.y = element_text(size = 30),
  legend.title = element_text(size = 30),
 # Adjust legend text size
) + theme(plot.caption = element_blank()) 
  theme(plot.caption = element_blank()) 
quartz()
total_phy
#pdf saved 13.38x8.24

############### ALL SOIL SAMPLES ###############################
##################################################################

# subset metadata to only include soil samples
soil = metadata[metadata$Samples %in% c('soil', 'control near', 'control far'),]
soilMeta = subset(metadata, Samples %in% c('soil', 'control near', 'control far'))

##### model selection for all diversity measures
#######################################################################
# use Samples column (differentiates between "soil", "control near" "control far")
#####  shannon ######
model_S1 = lm(shannon ~ Depth + Nest + Samples, data = soilMeta) 
check_model(model_S1)
summary(model_S1)
Anova(model_S1)

S2= lm(shannon ~ Nest + Depth * Samples, data = soilMeta) 
check_model(S2)
summary(S2)
Anova(S2)

S3 = lm(shannon ~ Depth * Nest * Samples, data = soilMeta)
check_model(S3)
summary(S3)
Anova(S3)

S4 = lm(shannon ~ Nest * Depth + Samples, data = soilMeta)
check_model(S4)
summary(S4)
Anova(S4)

S5 = lm(shannon ~ Nest * Samples + Depth, data = soilMeta)
check_model(S5)
summary(S5)
Anova(S5)

##### Faith's pd #####
model_F1 = lm(faith_pd ~ Depth + Nest + Samples, data = soilMeta) # ***
check_model(model_F1)
summary(model_F1)
Anova(model_F1)

F2= lm(faith_pd ~ Nest + Depth * Samples, data = soilMeta)
check_model(F2)
summary(F2)
Anova(F2)

F3 = lm(faith_pd ~ Depth * Nest * Samples, data = soilMeta)
check_model(F3)
summary(F3)
Anova(F3)

F4 = lm(faith_pd ~ Nest * Depth + Samples, data = soilMeta)
check_model(F4)
summary(F4)
Anova(F4)

F5 = lm(faith_pd ~ Nest * Samples + Depth, data = soilMeta)
check_model(F5)
summary(F5)
Anova(F5)

##### evenness #####
model_P1 = lm(pielou_evenness ~ Depth + Nest + Samples, data = soilMeta) # ***
check_model(model_P1)
summary(model_P1)
Anova(model_P1)

P2= lm(pielou_evenness ~ Nest + Depth * Samples, data = soilMeta)
check_model(P2)
summary(P2)
Anova(P2)

P3 = lm(pielou_evenness ~ Depth * Nest * Samples, data = soilMeta)
check_model(P3)
summary(P3)
Anova(P3)

P4 = lm(pielou_evenness ~ Nest * Depth + Samples, data = soilMeta)
check_model(P4)
summary(P4)
Anova(P4)

P5 = lm(pielou_evenness ~ Nest * Samples + Depth, data = soilMeta)
check_model(P5)
summary(P5)
Anova(P5)

##### observed ######
model_Obs1 = lm(observed_otus ~ Depth + Nest + Samples, data = soilMeta) # ***
check_model(model_Obs1)
summary(model_Obs1)
Anova(model_Obs1)

Obs2= lm(observed_otus ~ Nest + Depth * Samples, data = soilMeta)
check_model(Obs2)
summary(Obs2)
Anova(Obs2)

Obs3 = lm(observed_otus ~ Depth * Nest * Samples, data = soilMeta)
check_model(Obs3)
summary(Obs3)
Anova(Obs3)

Obs4 = lm(observed_otus ~ Nest * Depth + Samples, data = soilMeta)
check_model(Obs4)
summary(Obs4)
Anova(Obs4)

Obs5 = lm(observed_otus ~ Nest * Samples + Depth, data = soilMeta)
check_model(Obs5)
summary(Obs5)
Anova(Obs5)

# choose best fit models for each diversity index 
compare_performance(model_S1,S2,S3,S4,S5, rank = TRUE, verbose = FALSE) # model_S1
compare_performance(model_F1,F2,F3,F4,F5, rank = TRUE, verbose = FALSE) # model_F1
compare_performance(model_P1,P2,P3,P4,P5, rank = TRUE, verbose = FALSE) # model_P1
compare_performance(model_Obs1,Obs2,Obs3,Obs4,Obs5, rank = TRUE, verbose = FALSE) # model_obs1

### post hoc for all diversity measures
emmeans(model_S1, list(pairwise ~ Nest), adjust = "tukey") # shannon
emmeans(model_F1, list(pairwise ~ Nest), adjust = "tukey") # faith 
emmeans(model_P1, list(pairwise ~ Nest), adjust = "tukey") # evenness
emmeans(model_Obs1, list(pairwise ~ Nest), adjust = "tukey") # observed

## post hoc for evenness by soil type - anova was sig
emmeans(model_P1, list(pairwise ~ Samples), adjust = "tukey") # evenness 

############ PERMANOVA for ord plot #################################
####################################################################### 
# Extract sample data and distance matrix
physeq_soil_control = subset_samples(phy, Sample_Type %in% c("soil", "control"))
physeq_soil_control  =  physeq_soil_control %>%
  tax_fix() %>%
  phyloseq_validate()
physeq_soil_control  =  physeq_soil_control %>% tax_fix(unknowns = c("uncultured"))
physeq_soil_control@sam_data$Nest  =  factor(physeq_soil_control@sam_data$Nest, 
                                             levels = c("Apple", "Banana", "Carrot", "DC", "Eggplant"),
                                             labels = c("A", "B", "C", "D", "E")) 

sample_data2 = as(sample_data(physeq_soil_control), "data.frame")
bray_dist2 = distance(physeq_soil_control, method = "bray")

# Perform PERMANOVA
permanova_result_soilType = adonis2(bray_dist2 ~ Nest + Control_Sample1, data = sample_data2, permutations = 999)
print(permanova_result_soilType)

# perform pairwise permanova tests
pairwise_result_2a = pairwise.adonis(bray_dist2, sample_data2$Nest)
pairwise_result_2a

pairwise_result_2b = pairwise.adonis(bray_dist2, sample_data2$Control_Sample1)
pairwise_result_2b

############ FIGURE 5 - boxplot by colony and pcoa by colony ############
##########################################################################
#################### SHANNON PLOT #################
max_value = max(soilMeta$shannon) # Determine the maximum value of the y-axis to place the labels uniformly
label_y_position = max_value + 0.5
labels = data.frame(
  Nest = c("Apple", "Banana", "Carrot", "DC", "Eggplant"),
  shannon = rep(label_y_position, 5),  # Use the same y position for all labels 
  label = c("A", "B", "BC", "BC", "AC")
) 
ShannonFig5 <- ggboxplot(soilMeta, x = "Nest", y = "shannon", color = "Nest") + scale_x_discrete(labels = c("A", "B", "C", "D", "E")) + ggtitle("Shannon")
ShannonFig5 = ShannonFig5 + geom_boxplot(fill = c("red", "goldenrod1", "orange", "purple", "darkslateblue")) + theme_mine() + ggtitle("Shannon") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 35), legend.position = "none") + theme(axis.title.x = element_blank()) +theme(axis.text.x=element_text(size=30, face = "bold", angle = 0), axis.text.y = element_text(size = 25, angle = 90), axis.title.y = element_text(size = 30, angle = 90)) + 
  geom_text(data = labels, aes(x = Nest, y = shannon, label = label), 
            vjust = -0.5, size = 15) + ylim(NA, label_y_position + 0.5) + # set y-axis limits
  ylab("Shannon")
ShannonFig5

#################### FAITH PLOT #################

FaithFig5 <- ggboxplot(soilMeta, x = "Nest", y = "faith_pd", color = "Nest") + scale_x_discrete(labels = c("A", "B", "C", "D", "E")) + ggtitle("Faith PD")
FaithFig5 = FaithFig5 + geom_boxplot(fill = c("red", "goldenrod1", "orange", "purple", "darkslateblue")) + theme_mine() + ggtitle("Faith PD") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 35), legend.position = "none") + theme(axis.title.x = element_blank()) +theme(axis.text.x=element_text(size=30, face = "bold", angle = 0), axis.text.y = element_text(size = 25, angle = 90), axis.title.y = element_text(size = 30, angle = 90)) +  ylab("Faith's PD")

FaithFig5
                                                                                                                                       
#################### EVENNESS PLOT #################
max_value = max(soilMeta$pielou_evenness) # Determine the maximum value of the y-axis to place the labels uniformly
label_y_position = max_value + 0.01
labels = data.frame(
  Nest = c("Apple", "Banana", "Carrot", "DC", "Eggplant"),
  pielou_evenness = rep(label_y_position, 5),  # Use the same y position for all labels 
  label = c("A", "B", "B", "A", "A")
) 
EvenFig5 <- ggboxplot(soilMeta, x = "Nest", y = "pielou_evenness", color = "Nest") + scale_x_discrete(labels = c("A", "B", "C", "D", "E")) + ggtitle("Evenness")
EvenFig5 = EvenFig5 + geom_boxplot(fill = c("red", "goldenrod1", "orange", "purple", "darkslateblue")) + theme_mine() + ggtitle("Evenness") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 35), legend.position = "none") + theme(axis.title.x = element_blank()) + theme(axis.text.x=element_text(size=30, face = "bold", angle = 0), axis.text.y = element_text(size = 25, angle = 90), axis.title.y = element_text(size = 30, angle = 90)) + 
  geom_text(data = labels, aes(x = Nest, y = pielou_evenness, label = label), 
            vjust = -0.5, size = 15) + ylim(NA, label_y_position + 0.03) + # set y-axis limits
  ylab("Pielou's evenness")  
EvenFig5

#################### Observed PLOT #################
max_value = max(soilMeta$observed_otus) # Determine the maximum value of the y-axis to place the labels uniformly
label_y_position = max_value + 1
labels = data.frame(
  Nest = c("Apple", "Banana", "Carrot", "DC", "Eggplant"),
  observed_otus = rep(label_y_position, 5),  # Use the same y position for all labels 
  label = c("A", "B", "AB", "B", "AB")
) 
ObsFig5 <- ggboxplot(soilMeta, x = "Nest", y = "observed_otus", color = "Nest") + scale_x_discrete(labels = c("A", "B", "C", "D", "E")) + ggtitle("Observed ASVs")
ObsFig5 = ObsFig5 + geom_boxplot(fill = c("red", "goldenrod1", "orange", "purple", "darkslateblue")) + theme_mine() + ggtitle("Observed ASVs") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 35), legend.position = "none") + theme(axis.title.x = element_blank()) +theme(axis.text.x=element_text(size=30, face = "bold",angle = 0), axis.text.y = element_text(size = 25, angle = 90), axis.title.y = element_text(size = 30, angle = 90)) + 
  geom_text(data = labels, aes(x = Nest, y = observed_otus, label = label), 
            vjust = -0.5, size = 15) + ylim(NA, label_y_position + 35) + # set y-axis limits
  ylab("Observed")
ObsFig5

grid.arrange(ShannonFig5, FaithFig5, EvenFig5, ObsFig5, nrow= 2)
#pdf saved 16x11

############### FIGURE 5 PCOA by Colony ####################
##############################################################
## ord of total for mesa vs lowlands
#shape control near, control far, nest soil sample
ord_shape = physeq_soil_control %>%
  tax_transform(rank = "unique", trans = "identity") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Nest", fill = "Nest",
    shape = "Control_Sample1", alpha = 0.8,
    size = 2
  ) + 
  scale_shape_girafe_filled() +  ggplot2::stat_ellipse(
    ggplot2::aes(colour = Nest)
  ) + scale_colour_manual(values = c("red", "goldenrod1", "orange", "purple", "darkslateblue"))

ord_shape = ord_shape + scale_fill_manual(values = c("red", "goldenrod1", "orange", "purple", "darkslateblue")) + theme_mine() + theme(legend.position = c(0.95, 0.01),                                                                               legend.justification = c(0, 0)) + 
  theme(plot.caption = element_blank(),                                                                                       legend.title = element_text(size = 30),
          axis.text = element_text(size = 30),        
          legend.text = element_text(size = 25)                                                                     
        ) + 
    labs(shape = "Soil Type") # change legend title for shap
ord_shape <- ord_shape +
  theme(
    legend.position = "right",  # Move legend to the right
    legend.justification = "bottom"  # Center the legend
    )
# Mirror the plot by reversing the x-axis
ord_shape <- ord_shape + scale_y_reverse() 
ord_shape
#pdf saved 13.38x8.24

############# FIGURE 6 - EVENNESS PLOT for soil type ############# 
###################################################################
soilMeta$Samples = factor(soilMeta$Samples, 
                          levels = c("soil", "control near", "control far"),
                          labels = c("nest soil", "control near", "control far")) 

EvenFig6 <- ggboxplot(soilMeta, x = "Samples", y = "pielou_evenness", color = "Samples")  + ggtitle("Evenness")

max_value = max(soilMeta$pielou_evenness) # Determine the maximum value of the y-axis to place the labels uniformly
label_y_position = max_value + 0.01
labels = data.frame(
  Samples = c("nest soil", "control near", "control far"),
  pielou_evenness = rep(label_y_position, 3),  # Use the same y position for all labels 
  label = c("A", "AB", "B")
) 

EvenFig6L = EvenFig6 + geom_boxplot(fill = c("goldenrod3", "bisque3", "bisque1")) + theme_mine() + ggtitle("Evenness") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), legend.position = "none") + theme(axis.title.x = element_blank()) +theme(axis.text.x=element_text(size=30, face = "bold")) + 
  
  geom_text(data = labels, aes(x = Samples, y = pielou_evenness, label = label), 
            vjust = -0.5, size = 10) + ylim(NA, label_y_position + 0.015) + # set y-axis limits
  ylab("Pielou's evenness") + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y  = element_text(size = 30) )
  EvenFig6L

# control near is bisque3, control far is bisque1, soil samples is goldenrod3 

################## FIGURE 6B PCOA  BY SOIL TYPE ###############################
################################################################################ 
physeq_soil_control = subset_samples(phy, Sample_Type %in% c("soil", "control"))
physeq_soil_control <- physeq_soil_control %>%
  tax_fix() %>%
  phyloseq_validate()
physeq_soil_control <- physeq_soil_control %>% tax_fix(unknowns = c("uncultured"))
physeq_soil_control@sam_data$Control_Sample1 <- factor(physeq_soil_control@sam_data$Control_Sample1, 
                                            levels = c("nest_soil", "control near", "control far"),
                                            labels = c("nest soil", "control near", "control far"))
# ord plot by soil samples 
  all_soil <- physeq_soil_control %>%
    tax_transform(rank = "unique", trans = "identity") %>%
    dist_calc(dist = "bray") %>%
    ord_calc(
      method = "auto"
    ) %>% 
    ord_plot(
      axes = c(1, 2),
      colour = "Control_Sample1", fill = "Control_Sample1",
      shape = "circle", alpha = 0.8,
      size = 3
    ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(colour = Control_Sample1)
    )
  # add custon theme and color scales
  all_soil <- all_soil + theme_mine() + scale_colour_manual(values = c("goldenrod1", "bisque3", "bisque1")) + guides(colour = guide_legend(title = "Soil type"), fill = guide_legend(title = "Soil type")) + 
    theme(
      legend.position = "right", # adjust legend 
      legend.justification = c(0, 0), 
      legend.text = element_text(size = 25),
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      legend.title = element_text(size = 25),
    ) + theme(plot.caption = element_blank())
  all_soil = all_soil + scale_y_reverse()
  all_soil

############### NEST SOIL SAMPLES - no controls ############### #################################################################
### Run models for all diversity samples  
# subset data to only include nest soil samples
data <- read.csv("updated_clean_sample_metadata_diversity_metrics_239_1.csv")
data_soil = filter(data, Sample_Type == "soil")
data_soil1 = filter(data_soil, Chamber.Type.simple != "empty")

data_soil1$Chamber.Type.simple = factor(data_soil1$Chamber.Type.simple, 
                    levels =c ("ants", "brood", "brood/reproductive", "brood/seeds", "seeds"))
                   
###### shannon #####
model_lmS1 = lm(shannon ~ Nest + Depth + Chamber.Type.simple, data = data_soil1)
check_model(model_lmS1)
summary(model_lmS1)
Anova(model_lmS1)

lmS2 = lm(shannon ~ Nest * Depth * Chamber.Type.simple, data = data_soil1) 
check_model(lmS2)
summary(lmS2)
Anova(lmS2)

lmS3 = lm(shannon ~ Nest * Depth + Chamber.Type.simple, data = data_soil1) 
check_model(lmS3)
summary(lmS3)
Anova(lmS3)

lmS4 = lm(shannon ~ Nest + Depth * Chamber.Type.simple, data = data_soil1) 
check_model(lmS4)
summary(lmS4)
Anova(lmS4)

lmS5 = lm(shannon ~ Depth + Nest * Chamber.Type.simple, data = data_soil1) 
check_model(lmS5)
summary(lmS5)
Anova(lmS5)

##### faith #####
model_lmF1 = lm(faith_pd ~ Nest + Depth + Chamber.Type.simple, data = data_soil1)
check_model(model_lmF1)
summary(model_lmF1)
Anova(model_lmF1)

lmF2 = glm(faith_pd ~ Nest * Depth * Chamber.Type.simple, data = data_soil1) #****
check_model(lmF2)
summary(lmF2)
Anova(lmF2)

lmF3 = lm(faith_pd ~ Nest * Depth + Chamber.Type.simple, data = data_soil1)
check_model(lmF3)
summary(lmF3)
Anova(lmF3)

lmF4 = lm(faith_pd ~ Nest + Depth * Chamber.Type.simple, data = data_soil1)
check_model(lmF4)
summary(lmF4)
Anova(lmF4)

lmF5 = lm(faith_pd ~ Depth + Nest * Chamber.Type.simple, data = data_soil1) 
check_model(lmF5)
summary(lmF5)
Anova(lmF5)

##### evenness #####
model_lmP1 = lm(pielou_evenness ~ Nest + Depth + Chamber.Type.simple, data = data_soil1)
check_model(model_lmP1)
summary(model_lmP1)
Anova(model_lmP1) 

lmP2 = lm(pielou_evenness ~ Nest * Depth * Chamber.Type.simple, data = data_soil1) 
check_model(lmP2)
summary(lmP2)
Anova(lmP2)

lmP3 = lm(pielou_evenness ~ Nest * Depth + Chamber.Type.simple, data = data_soil1)
check_model(lmP3)
summary(lmP3)
Anova(lmP3)

lmP4 = lm(pielou_evenness ~ Nest + Depth * Chamber.Type.simple, data = data_soil1)
check_model(lmP4)
summary(lmP4)
Anova(lmP4)

lmP5 = lm(pielou_evenness ~ Depth + Nest * Chamber.Type.simple, data = data_soil1) 
check_model(lmP5)
summary(lmP5)
Anova(lmP5)

##### observed #####
model_lmO1 = lm(observed_otus ~ Nest + Depth + Chamber.Type.simple, data = data_soil1)
check_model(model_lmO1)
summary(model_lmO1)
Anova(model_lmO1) 

lmO2 = lm(observed_otus ~ Nest * Depth * Chamber.Type.simple, data = data_soil1)
check_model(lmO2)
summary(lmO2)
Anova(lmO2)

lmO3 = lm(observed_otus ~ Nest * Depth + Chamber.Type.simple, data = data_soil1)
check_model(lmO3)
summary(lmO3)
Anova(lmO3)

lmO4= lm(observed_otus ~ Nest + Depth * Chamber.Type.simple, data = data_soil1)
check_model(lmO4)
summary(lmO4)
Anova(lmO4)

lmO5 = lm(observed_otus ~ Depth + Nest * Chamber.Type.simple, data = data_soil1) 
check_model(lmO5)
summary(lmO5)
Anova(lmO5)

#find best fit model
compare_performance(model_lmS1, lmS2, lmS3, lmS4, lmS5, rank = TRUE, verbose = FALSE) # model_lmS1
compare_performance(model_lmF1, lmF2, lmF3, lmF4,lmF5,  rank = TRUE, verbose = FALSE) # model_lmF1
compare_performance(model_lmP1, lmP2, lmP3, lmP4, lmP5, rank = TRUE, verbose = FALSE) # model_lmP1
compare_performance(model_lmO1, lmO2, lmO3, lmO4,lmO5, rank = TRUE, verbose = FALSE) # model_lmO1

##### post hoc for all diversity measures #####

emmeans(model_lmS1, list(pairwise ~ Chamber.Type.simple), adjust = "tukey") # shannon
emmeans(model_lmF1, list(pairwise ~ Chamber.Type.simple), adjust = "tukey") # faith 
emmeans(model_lmP1, list(pairwise ~ Chamber.Type.simple), adjust = "tukey") # evenness
emmeans(model_lmO1, list(pairwise ~ Chamber.Type.simple), adjust = "tukey") # observed


########## PERMANOVA #####################
############################################ 
# extract sample data and distance matrix
sample_data3 = as(sample_data(physeq_soil_NE), "data.frame")
bray_dist3 = distance(physeq_soil_NE, method = "bray")
# perform PERMANOVA
permanova_result_nestSoil = adonis2(bray_dist3~Nest+Chamber.Type.simple,data = sample_data3, permutations = 999 )
print(permanova_result_nestSoil)

# perform pairwise permanova tests
pairwise_result3_nest <- pairwise.adonis(bray_dist3, sample_data3$Nest)
pairwise_result3_nest

######################## FIGURE 7A - BOX PLOTS BY CHAMBER ######################## 
###################################################################################   
data <- read.csv("updated_clean_sample_metadata_diversity_metrics_239_1.csv")
data_soil = filter(data, Sample_Type == "soil")
data_soil1 = filter(data_soil, Chamber.Type.simple != "empty")

data_soil1$Chamber.Type.simple = factor(data_soil1$Chamber.Type.simple, 
                                        levels =c ("ants", "brood", "brood/reproductive", "brood/seeds", "seeds"),
                                        labels = c("ants","brood", "br + rep", "br + seeds", "seeds"))

Shannon <- ggboxplot(data_soil1, x = "Chamber.Type.simple", y = "shannon", color = "Chamber.Type.simple") + ggtitle("Shannon")

Shannon = Shannon + geom_boxplot(fill = c("darkviolet", "mediumorchid1", "plum1", "cyan2", "olivedrab3")) + theme_mine() + ggtitle("Shannon") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), legend.position = "none") + theme(axis.title.x = element_blank()) +theme(axis.text.x=element_text(size=15, angle = 0), axis.text.y = element_text(size = 20, angle = 90), axis.title.y = element_text(size = 30, angle = 90)) +  ylab("Shannon") + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y  = element_text(size = 30) )
Shannon

faith_pd <- ggboxplot(data_soil1, x = "Chamber.Type.simple", y = "faith_pd", color = "Chamber.Type.simple") + ggtitle("Faith PD")
faith_pd = faith_pd + geom_boxplot(fill = c("darkviolet", "mediumorchid1", "plum1", "cyan2", "olivedrab3")) + theme_mine() + ggtitle("Faith PD") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), legend.position = "none") + theme(axis.title.x = element_blank()) +theme(axis.text.x=element_text(size=15, angle = 0), axis.text.y = element_text(size = 20, angle = 90), axis.title.y = element_text(size = 30, angle = 90)) +  ylab("Faith's PD") + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y  = element_text(size = 30) )
faith_pd

Observed <- ggboxplot(data_soil1, x = "Chamber.Type.simple", y = "observed_otus", color = "Chamber.Type.simple") + ggtitle("Observed ASVs")
Observed = Observed + geom_boxplot(fill = c("darkviolet", "mediumorchid1", "plum1", "cyan2", "olivedrab3")) + theme_mine() + ggtitle("Observed ASVs") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), legend.position = "none") + theme(axis.title.x = element_blank()) +theme(axis.text.x=element_text(size=15, angle = 0), axis.text.y = element_text(size = 20, angle = 90), axis.title.y = element_text(size = 30, angle = 90)) +  ylab("Observed") + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y  = element_text(size = 30) )
Observed

Even <- ggboxplot(data_soil1, x = "Chamber.Type.simple", y = "pielou_evenness", color = "Chamber.Type.simple") + ggtitle("Evenness")
Even = Even + geom_boxplot(fill = c("darkviolet", "mediumorchid1", "plum1", "cyan2", "olivedrab3")) + theme_mine() + ggtitle("Evenness") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), legend.position = "none") + theme(axis.title.x = element_blank()) + theme(axis.text.x=element_text(size=18, angle = 0), axis.text.y = element_text(size = 20, angle = 90), axis.title.y = element_text(size = 30, angle = 90)) +  ylab("Pielou's evenness") + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y  = element_text(size = 30) )
Even

grid.arrange(Shannon, faith_pd, Even, Observed, nrow= 2)
#pdf saved 16x11


######## FIGURE 7 - bray plot for only in-nest soil by colored chamber type ########
#####################################################################################  
physeq_soil = subset_samples(phy, Sample_Type == "soil")
physeq_soil_NE = subset_samples(physeq_soil, Chamber.Type.simple != "empty")
physeq_soil_NE@sam_data$Chamber.Type.simple <- factor(physeq_soil_NE@sam_data$Chamber.Type.simple,
                          levels = c("ants", "brood", "brood/reproductive", "brood/seeds", "seeds"),
                          labels = c("ants", "brood", "br + rep", "br + seeds", "seeds"))

physeq_soil_NE@sam_data$Nest <- factor(physeq_soil_NE@sam_data$Nest, 
                      levels = c("Apple", "Banana", "Carrot", "DC", "Eggplant"),
                    labels = c("A", "B", "C", "D", "E"))


soil <-physeq_soil_NE %>%
  tax_transform(rank = "unique", trans = "identity") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Chamber.Type.simple", fill = "Chamber.Type.simple", 
    shape = "Nest",
    alpha = 0.8,
    size = 3
  ) +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = Chamber.Type.simple)
  )
soil1 <- soil + scale_colour_manual(values = c( "darkviolet", "mediumorchid1", "plum1", "cyan2", "olivedrab3")) + theme_mine() + guides(colour = guide_legend(title = "Chamber type"), fill = guide_legend(title = "Chamber type"))+
  theme(
    legend.position = "right", # adjust legend 
    legend.justification = c(0, 0), 
    legend.text = element_text(size = 25),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.title = element_text(size = 25),
    ) 
soil1 = soil1 + scale_y_reverse() + theme(plot.caption = element_blank())
soil1

### only chamber type 
soilCT <-physeq_soil_NE %>%
  tax_transform(rank = "unique", trans = "identity") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Chamber.Type.simple", fill = "Chamber.Type.simple",
    alpha = 0.8,
    size = 3
  ) +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = Chamber.Type.simple)
  )
soilCT2 <- soilCT + scale_colour_manual(values = c( "darkviolet", "mediumorchid1", "plum1", "cyan2", "olivedrab3")) + theme_mine() + guides(colour = guide_legend(title = "Chamber type"), fill = guide_legend(title = "Chamber type")) +
theme(
  legend.position = "right", # adjust legend 
  legend.justification = c(0, 0), 
  legend.text = element_text(size = 25),
  axis.title.x = element_text(size = 30),
  axis.title.y = element_text(size = 30),
  legend.title = element_text(size = 25),
) + scale_y_reverse()
