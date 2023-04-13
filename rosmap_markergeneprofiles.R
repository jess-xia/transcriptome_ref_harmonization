library(readr)
library(tibble)
library(knitr)
library(matrixStats)
#devtools::install_github('oganm/markerGeneProfile', force = T) # install marker gene profile tool from github
library(markerGeneProfile)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(cowplot)

# Load data

# Load metadata
meta <- read_csv("/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/raw_gene_counts/metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv") %>%
  # filter duplicated genes 
  filter(duplicated(specimenID) == FALSE) %>%
  dplyr::select(specimenID,
         age_death, 
         msex, 
         pmi, 
         rin = RIN)

# Combined batches to make large count matrix
counts_batch1 <- read_tsv("/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/raw_gene_counts/gene_quantification/rosmap_batch1_stranded/ROSMAP_batch1_gene_all_counts_matrix_clean.txt")
counts_batch2 <- read_tsv("/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/raw_gene_counts/gene_quantification/rosmap_batch2_stranded/ROSMAP_batch2_gene_all_counts_matrix_clean.txt")
counts <- cbind(counts_batch1, counts_batch2)

counts <- counts[5:nrow(counts), ] # Keep only gene counts
counts$feature <- gsub("\\..*", "", counts$feature) # remove periods and numbers in gene names used to avoid duplicates 
counts <- counts[duplicated(counts$feature) == FALSE,] # filter duplicated genes 
rownames(counts) <- counts$feature
counts <- subset(counts, select = -c(feature) )


# Normalize counts & process

# Convert counts to matrix, normalize using CPM, add 0.1 to every observation, and log2 transform 
counts_cpm = edgeR::cpm(data.matrix(counts), log = TRUE, prior.count = 0.1)
kable(counts_cpm[1:5, 1:10]) # preview
# calculate standard deviation per row
gene_sds = rowSds(counts_cpm, na.rm = T) # missing values are excluded
gene_mat = counts_cpm[gene_sds > .1, ] %>% # keep rows (genes) with sd greater than 0.1 across all samples 
  as.data.frame() %>%
  rownames_to_column(var = "gene_symbol")



# Cell-proportion estimation
marker_data <- readRDS("/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/M1_all_group_markers.rds")
cell_types = marker_data$cluster %>% unique()

# organize markers into a list, which is the format will need for next steps 
marker_list = lapply(cell_types, function(cell_type){
  return(marker_data %>% filter(cluster == cell_type) %>% pull(gene) %>% unlist())
})
names(marker_list) = cell_types # Set names as the cell types

# Run marker gene profile (MGP) analysis
estimations =  mgpEstimate(
  exprData = gene_mat,
  genes = marker_list,
  geneColName = 'gene_symbol',
  outlierSampleRemove = FALSE, # should outlier samples removed. This is done using boxplot stats
  geneTransform = NULL, # this is the default option for geneTransform
  groups = NULL, # if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority = TRUE)


# get proportion estimates as data frame
estimations_df = as.data.frame(estimations$estimates) %>%
  rownames_to_column(var = "specimenID")
saveRDS(estimations_df, "/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/rosmap_mgp_estimations.rds")

# proportions are unit-less and negative numbers are just a result of PCA
# we can scale proportions between 0-1 for visualization purposes 
scale0 = function(x){
  (x-min(x))/(max(x)-min(x))
}
estimations_df_norm = scale0(estimations_df[,-1]) # don't include the first column
estimations_df_norm = estimations_df_norm %>%
  mutate(specimenID = estimations_df$specimenID) %>%
  dplyr::select(specimenID, everything())


# merge cell type proportions with sample metadata
mgp_df = inner_join(meta, estimations_df_norm, by = "specimenID") %>%
  # pivot longer so there is one row per sample and cell type
  # this type of data structure is preferred by ggplot
  pivot_longer(-colnames(meta),
               names_to = "cell_type",
               values_to = "cell_proportion")
# fix labels 
mgp_df$cell_type = gsub("\\.", "/", mgp_df$cell_type)
# preview 
kable(mgp_df[1:10,])


# Quality-control of cell-proportion estimates
# loop through each cell type 
for(i in 1:length(cell_types)){
  # get the expression values of markers kept by mgpEstimate() per bulk sample 
  cells_df = estimations$usedMarkerExpression[i] %>% as.data.frame()
  # get list of markers kept by mgpEstimate()
  masterlist = paste0(rownames(cells_df), collapse=', ')
  # number of markers kept by mgpEstimate()
  num_markers = length(rownames(cells_df))
  # ratio of markers removed
  rm_marker_ratios = estimations$removedMarkerRatios[i]
  # check if the trimmed PCs are not(!) empty (NULL)
  if(!is.null(estimations$trimmedPCAs[[i]])){
    # get the percent variance that each principal component (PC) captures in the data 
    percent_variance = ((summary(estimations$trimmedPCAs[[i]]))[6]) %>% as.data.frame()
    # get the percent variance that the first principal component (PC1) captures in the data 
    percent_variance_PC1 = percent_variance[2,1]
  }
  else{
    # otherwise, the set this to NA
    percent_variance_PC1 = NA
  }
  # build a dataframe with all our metrics
  # initialize dataframe if this is the first iteration of the for() loop
  if(i==1){
    master_df = data.frame("markers_used" = masterlist, 
                           "removed_marker_ratios" = rm_marker_ratios,
                           "percent_variance_PC1" = percent_variance_PC1, 
                           "num_markers" = num_markers)  
  }
  # bind previous iterations results
  else{
    df = data.frame("markers_used" = masterlist, 
                    "removed_marker_ratios" = rm_marker_ratios,
                    "percent_variance_PC1" = percent_variance_PC1, 
                    "num_markers" = num_markers)
    master_df = rbind(master_df, df)
  }
}
QC_metrics = rownames_to_column(master_df, var = "celltype")

# Plot the number of markers per cell type
QC_metrics %>%  ggplot(aes(x = celltype, y = num_markers)) +
  theme_minimal() +
  geom_bar(stat = "identity", fill = "#e0abf5") +
  geom_hline(yintercept = 4) + 
  labs(title="Plot of Number of Markers Per Celltype", 
       x="Cell Type", y = "Markers Used")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  coord_flip()




# Plot the percent variance explained by PC (principal component) 1
QC_metrics %>%  ggplot(aes(x=celltype, y=percent_variance_PC1))+
  geom_bar(stat = "identity", fill =ifelse(QC_metrics$percent_variance_PC1 > 0.35, "#AFEEEE", "#808080")) +
  geom_hline(yintercept = 0.35) +
  labs(title="Percent Variance Explained by Each MGP",
       x="MGPs", y = "Percent Variance Explained by PC 1")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90)) + coord_flip()


# Run MGP with markers from Micaela's paper
pub_markers = read.csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'))[,-1]
pub_cell_types = pub_markers$Subclass %>% unique()

# organize markers into a list, which is the format will need for next steps 
pub_marker_list = lapply(pub_cell_types, function(cell_type){
  return(pub_markers %>% filter(Subclass == cell_type) %>% pull(Ensembl.gene.ID) %>% unlist())
})
names(pub_marker_list) = pub_cell_types # Set names as the cell types
pub_marker_list <- strsplit(pub_marker_list, ", ")
# Split ENSG gene codes that are separated by commas
pub_marker_list_reformatted <- lapply(pub_marker_list, function(sub_list) {
  unlist(lapply(sub_list, function(string) {
    strsplit(string, ", ")
  }))
})
# Remove NA values from all sublists using lapply
pub_marker_list_reformatted <- lapply(pub_marker_list_reformatted, function(sub_list) {
  unlist(lapply(sub_list, function(x) na.omit(x)))
})

# Run marker gene profile (MGP) analysis
pub_estimations =  mgpEstimate(
  exprData = gene_mat,
  genes = pub_marker_list,
  geneColName = 'gene_symbol',
  outlierSampleRemove = FALSE, # should outlier samples removed. This is done using boxplot stats
  geneTransform = NULL, # this is the default option for geneTransform
  groups = NULL, # if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority = TRUE)


# get proportion estimates as data frame
pub_estimations_df = as.data.frame(pub_estimations$estimates) %>%
  rownames_to_column(var = "specimenID")
saveRDS(pub_estimations_df, "/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/pub_rosmap_mgp_estimations.rds")



