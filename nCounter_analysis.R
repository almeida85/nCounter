#
# This code analyzes the gene expression data generated with the NanoString nCounter assay 
# for five patients in two timepoints.
# First, the data is normalized and a heatmap is generated for the positive and negtive control genes.
# Second, a summary statistics boxplot is generated.
#

# Loading necessary packages and libraries
require('NanoStringNorm')
library(NanoStringNorm)
library(gplots)
library(RColorBrewer)
library(scales)
library(stats)
library(ggplot2)
library(reshape2)

# FUNCTIONS

norm <- function(raw_data) {
  # Data normalization
  data_norm <- NanoStringNorm(raw_data,
                              CodeCount = 'geo.mean',
                              Background = 'mean.2sd',
                              SampleContent = 'housekeeping.geo.mean',
                              round.values = TRUE,
                              take.log = TRUE)
  
  return(data_norm)
}

get_graphs <- function(genes_norm, gene_list){
  # This function takes the normalized genes and a list of gene to study
  # and return an boxplot of those genes
  
  if (nargs() < 2 | length(gene_list) == 0) {
    cat("ERROR: Gene lists is empty or missing!!!")
    return(0)
  }
  
  genes_to_plot <- list()
  
  # Time points...
  treatments <- c("Baseline", "Post_treatment")
  
  # Getting the genes to plot
  for (i in seq_along(gene_list)) {
    genes_to_plot[[i]] <- subset(genes_norm, Name == gene_list[i])
  }
  
  # Treatment groups for each gene
  geneA_baseline <- c(genes_to_plot[[1]]$geneA_patient1_baseline_mRNA,
                     genes_to_plot[[1]]$geneA_patient2_baseline_mRNA,
                     genes_to_plot[[1]]$geneA_patient3_baseline_mRNA,
                     genes_to_plot[[1]]$geneA_patient4_baseline_mRNA,
                     genes_to_plot[[1]]$geneA_patient5_baseline_mRNA)

  geneA_post_treatment <- c(genes_to_plot[[1]]$geneA_patient1_post_mRNA,
                           genes_to_plot[[1]]$geneA_patient2_post_mRNA,
                           genes_to_plot[[1]]$geneA_patient3_post_mRNA,
                           genes_to_plot[[1]]$geneA_patient4_post_mRNA,
                           genes_to_plot[[1]]$geneA_patient5_post_mRNA)
  
  geneB_baseline <- c(genes_to_plot[[2]]$geneB_patient1_baseline_mRNA,
                     genes_to_plot[[2]]$geneB_patient2_baseline_mRNA,
                     genes_to_plot[[2]]$geneB_patient3_baseline_mRNA,
                     genes_to_plot[[2]]$geneB_patient4_baseline_mRNA,
                     genes_to_plot[[2]]$geneB_patient5_baseline_mRNA)
  
  geneB_post_treatment <- c(genes_to_plot[[2]]$geneB_patient1_post_mRNA,
                           genes_to_plot[[2]]$geneB_patient2_post_mRNA,
                           genes_to_plot[[2]]$geneB_patient3_post_mRNA,
                           genes_to_plot[[2]]$geneB_patient4_post_mRNA,
                           genes_to_plot[[2]]$geneB_patient5_post_mRNA)
  
  geneA_timepoints_groups <- data.frame(geneA_baseline, geneA_post_treatment)
  names(geneA_timepoints_groups) <- treatments
  geneA_timepoints_groups$Gene = gene_list[1]
  
  geneB_timepoints_groups <- data.frame(geneB_baseline, geneB_post_treatment)
  names(geneB_timepoints_groups) <- treatments
  geneB_timepoints_groups$Gene = gene_list[2]
  
  # Merging the gene tables and adjusting column names
  geneA_geneB_table <- melt(rbind(geneA_timepoints_groups, geneB_timepoints_groups),
                           id.vars='Gene', 
                           measure.vars = treatments)
  colnames(geneA_geneB_table)[2] <- "Treatment"
  colnames(geneA_geneB_table)[3] <- "Norm_value"
  
  # Plotting
  ggplot(geneA_geneB_table, aes(x=Gene, y=Norm_value, fill = Treatment)) +

    stat_summary(fun.y = mean,
                 geom = "point",
                 size = 5,
                 shape = 25,
                 color = "black",
                 position = position_dodge(width=0.7)) +

    geom_boxplot(outlier.color = "blue",
                 outlier.shape = 4,
                 outlier.size = 3,
                 alpha = 0.5) + # <- Transparency

    geom_jitter(position = position_jitterdodge(jitter.width = .1, dodge.width = 1),
                size = 3) +
    theme_bw()
  
  return(genes_to_plot)
}


# Read RCC file(s)
RCC_files_raw <- read.markup.RCC()

RCC_norm <- norm(RCC_files_raw)

genes_norm <- RCC_norm$normalized.data

# Write the list of the genes to study
gene_list <- list("geneA", "geneB")

get_graphs(genes_norm, gene_list)

# Select the Table with genes information
#gene_tables <- RCC_files_raw$x




