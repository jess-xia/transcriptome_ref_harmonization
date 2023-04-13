# Comparing mgp predicted proportions with ROSMAP sn. 
# Both the MGP I generated and from Micaela's paper

library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

# Import sn proportions and meta data
sn_proportions_raw <- read_csv("/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/rosmap_single_nuc_proportions.csv")
meta_raw <- read_csv("/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/raw_gene_counts/metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv") 

# Remove duplicates for mapping specimenID to projid
matcher <- meta_raw[duplicated(meta_raw$specimenID) == FALSE, c("specimenID", "projid")]
matcher <- matcher[duplicated(matcher$projid) == FALSE,]

# Assign specimenID to sn proportions which only has projid
sn_proportions <- merge(sn_proportions_raw, matcher, by.x = "ID", by.y = "projid", all.x = TRUE, all.y = FALSE)
sn_proportions <- sn_proportions[!is.na(sn_proportions$specimenID),] %>%
  subset(select = -c(X1))
colnames(sn_proportions) <- make.names(colnames(sn_proportions))

# Prepare data individually
rosmap_mgp_estimations <- readRDS("/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/rosmap_mgp_estimations.rds")
mgp_sst <- rosmap_mgp_estimations %>%
  subset(select = c(specimenID, Sst))
colnames(mgp_sst)[2] <- "mgp_sst"
sn_sst <- sn_proportions[c("specimenID", "Sst")]
colnames(sn_sst)[2] <- "sn_sst"

# Combined data
combined_sst <- merge(sn_sst, mgp_sst, by = "specimenID")
combined_sst <- combined_sst[, -1]

# Plot
ggplot(combined_sst, aes(x=sn_sst, y=mgp_sst)) +
  geom_point(size=2, shape=23) +
  geom_smooth(method = "lm", se = F) +
  stat_cor(method = "pearson")


# Match up column names for sn_proportions & rosmap_mgp_estimations
setdiff(colnames(sn_proportions), colnames(rosmap_mgp_estimations))
setdiff(colnames(rosmap_mgp_estimations), colnames(sn_proportions))
common_subclass <- intersect(colnames(sn_proportions), colnames(rosmap_mgp_estimations))
common_subclass <- common_subclass[common_subclass != "specimenID"]

rosmap_mgp_estimations <- dplyr::rename(rosmap_mgp_estimations, Micro = Micro.PVM)
# sn_proportions <- dplyr::rename()

# Create an empty list to store the plots
plot_list <- list()

for (i in seq(length(common_subclass))) {
  mgp <- rosmap_mgp_estimations[,c("specimenID", common_subclass[i])]
  colnames(mgp)[2] <- "mgp"
  sn <- sn_proportions[c("specimenID", common_subclass[i])]
  colnames(sn)[2] <- "sn"
  
  # Combined data
  combined <- merge(sn, mgp, by = "specimenID")
  combined <- combined[, -1]
  plot <- ggplot(combined, aes(x=sn, y=mgp)) +
    geom_point(size=0.5, shape=1) +
    geom_smooth(method = "lm", se = F) +
    stat_cor(method = "pearson") +
    ggtitle(common_subclass[i])
  
  # Add the plot to the list
  plot_list[[i]] <- plot
}

# Arrange the plots in a 2x2 grid
grid.arrange(grobs = plot_list)
  




# Plot mgp results generated from markers from Micaela's paper
pub_rosmap_mgp_estimations <- readRDS("/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/pub_rosmap_mgp_estimations.rds") %>%
  dplyr::rename(Astro = Astrocyte, 
                Endo = Endothelial, 
                Pvalb = PVALB, 
                Vip = VIP, 
                Sst = SST, 
                Pax6 = PAX6, 
                Lamp5 = LAMP5,
                Micro = Microglia, 
                Oligo = Oligodendrocyte, 
                L6.IT.Car3 = L5.6.IT.Car3)

sn_proportions$IT <- apply(sn_proportions[, c("L5.IT", "L6.IT", "L2.3.IT")], 1, sum)
sn_proportions <- sn_proportions[, !(names(sn_proportions) %in% c("L5.IT", "L6.IT", "L2.3.IT"))]

# Match up column names for sn_proportions & rosmap_mgp_estimations
setdiff(colnames(sn_proportions), colnames(pub_rosmap_mgp_estimations))
setdiff(colnames(pub_rosmap_mgp_estimations), colnames(sn_proportions))
common_subclass <- intersect(colnames(sn_proportions), colnames(pub_rosmap_mgp_estimations))
common_subclass <- common_subclass[common_subclass != "specimenID"]


plot_correlation(sn_proportions, pub_rosmap_mgp_estimations, "pub_mgp")









