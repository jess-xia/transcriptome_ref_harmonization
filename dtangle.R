# install.packages("dtangle")
library(dtangle)

# dtangle requires each row contains expression measurements for a particular sample. Each columm
# contains the measurements of the same gene over all individuals. 
# Combined batches to make large count matrix
counts_batch1 <- read_tsv("/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/raw_gene_counts/gene_quantification/rosmap_batch1_stranded/ROSMAP_batch1_gene_all_counts_matrix_clean.txt")
counts_batch2 <- read_tsv("/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/raw_gene_counts/gene_quantification/rosmap_batch2_stranded/ROSMAP_batch2_gene_all_counts_matrix_clean.txt")
counts <- cbind(counts_batch1, counts_batch2)

counts <- counts[5:nrow(counts), ] # Keep only gene counts
counts$feature <- gsub("\\..*", "", counts$feature) # remove periods and numbers in gene names used to avoid duplicates 
counts <- counts[duplicated(counts$feature) == FALSE,] # filter duplicated genes 
rownames(counts) <- counts$feature
counts <- subset(counts, select = -c(feature) )


# Convert counts to matrix, normalize using CPM, add 0.1 to every observation, and log2 transform 
counts_cpm = edgeR::cpm(data.matrix(counts), log = TRUE, prior.count = 0.1)
bulk_counts <- t(counts_cpm)

sc_allcounts_raw <- readRDS("~/transcriptome-ref-harmonization/sc_allcounts_raw.rds")
sc_allcounts <- column_to_rownames(sc_allcounts_raw, "Group.1") # Set first column as rownames, needed to convert to matrix
sc_allcounts <- sc_allcounts[,new_meta$cell] #filter for the cells we specified previously (removing outliers). may display error because cells in metadata are not in count matrix
new_meta <- new_meta[colnames(sc_allcounts),]
#new_meta <- new_meta[,-1] # Remove first column X in dataframe

# Convert dataframe to matrix 
new_count_matrix <- as.matrix(sc_allcounts, sparse = TRUE) 
# Convert counts to matrix, normalize using CPM, add 0.1 to every observation, and log2 transform 
new_counts_cpm = edgeR::cpm(new_count_matrix, log = TRUE, prior.count = 0.1)
kable(new_counts_cpm[1:5, 1:10]) # preview

azimuth <- read_tsv("/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/azimuth_pred.tsv") %>%
  dplyr::rename(azimuth_subclass = predicted.subclass) %>%
  subset(select = c(cell, azimuth_subclass)) %>%
  column_to_rownames("cell")
azimuth$azimuth_subclass <- make.names(azimuth$azimuth_subclass)


# t_counts_sc is the reference counts
t_counts_sc <- t(new_counts_cpm)
sc_ref <- as.data.frame(t_counts_sc) %>%
  merge(azimuth, ., by = "row.names") %>%
  tibble::rowid_to_column("index")
sc_ref <- sc_ref[, 1:3]

pure_samps <- list()
subclasses <- make.names(unique(sc_ref$azimuth_subclass))
for (i in seq(length(subclasses))){
  temp <- sc_ref$index[sc_ref$azimuth_subclass == subclasses[i]]
  pure_samps <- c(pure_samps, list(temp))
}
names(pure_samps) <- subclasses

# Make sure column names are the same
common_cols <- intersect(colnames(t_counts), colnames(t_counts_sc))
t_counts <- t_counts[, common_cols]
t_counts_sc <- t_counts_sc[, common_cols]



dt_out <- dtangle(Y = t_counts, reference = t_counts_sc, pure_samples = pure_samps)
# saveRDS(dt_out$estimates, "/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/rosmap_dt_estimations.rds")
rosmap_dt_estimations <- readRDS("/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/rosmap_dt_estimations.rds") %>%
  as.data.frame() %>%
  dplyr::rename(Micro = Micro.PVM) 

rosmap_dt_estimations <- rownames_to_column(rosmap_dt_estimations, "specimenID")


sn_proportions_raw <- read_csv("/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/rosmap_single_nuc_proportions.csv")

# Assign specimenID to sn proportions which only has projid
sn_proportions <- merge(sn_proportions_raw, matcher, by.x = "ID", by.y = "projid", all.x = TRUE, all.y = FALSE)
sn_proportions <- sn_proportions[!is.na(sn_proportions$specimenID),] %>%
  subset(select = -c(X1))
colnames(sn_proportions) <- make.names(colnames(sn_proportions))


# Match up column names for sn_proportions & rosmap_mgp_estimations
setdiff(colnames(sn_proportions), colnames(rosmap_dt_estimations))
setdiff(colnames(rosmap_dt_estimations), colnames(sn_proportions))
common_subclass <- intersect(colnames(sn_proportions), colnames(rosmap_mgp_estimations))
common_subclass <- common_subclass[common_subclass != "specimenID"]

plot_correlation <- function(sn_proportions, rosmap_estimation, deconv_type){
  # Create an empty list to store the plots
  plot_list <- list()
  
  for (i in seq(length(common_subclass))) {
    deconv <- rosmap_estimation[,c("specimenID", common_subclass[i])]
    colnames(deconv)[2] <- deconv_type
    sn <- sn_proportions[c("specimenID", common_subclass[i])]
    colnames(sn)[2] <- "sn"
    
    # Combined data
    combined <- merge(sn, deconv, by = "specimenID")
    combined <- combined[, -1]
    plot <- ggplot(combined, aes(x=sn, y=get(deconv_type))) +
      geom_point(size=0.5, shape=1) +
      geom_smooth(method = "lm", se = F) +
      stat_cor(method = "pearson") +
      ggtitle(common_subclass[i]) + 
      labs(y = deconv_type)
    
    # Add the plot to the list
    plot_list[[i]] <- plot
  }
  
  # Arrange the plots in a 2x2 grid
  grid.arrange(grobs = plot_list)

}

plot_correlation(sn_proportions, rosmap_dt_estimations, "dt")
