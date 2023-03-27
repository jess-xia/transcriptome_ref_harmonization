library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

# Import sn proportions and meta data
sn_proportions_raw <- read_csv("/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/rosmap_single_nuc_proportions.csv")
meta_raw <- read_csv("/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/raw_gene_counts/metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv") 

# Remove duplicates for mapping specimenID to projid
matcher <- meta_raw[duplicated(meta_raw$specimenID) == FALSE, c("specimenID", "projid")]
matcher <- matcher[duplicated(matcher$projid) == FALSE,]

# Assign specimenID to sn proportions which only has projid
sn_proportions <- merge(sn_proportions_raw, matcher, by.x = "ID", by.y = "projid", all.x = TRUE, all.y = FALSE)
sn_proportions <- sn_proportions[!is.na(sn_proportions$specimenID),]

# Prepare data individually
rosmap_mgp_estimations <- readRDS("/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/rosmap_mgp_estimations.rds")
mgp_sst <- rosmap_mgp_estimations[c("specimenID", "SST")]
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

