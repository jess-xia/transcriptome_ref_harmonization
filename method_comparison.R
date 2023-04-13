library(tidyr)
library(ggplot2)

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


# Plot mgp results generated from markers from Micaela's paper
sn_proportions_forpub <- sn_proportions
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
# Set specimenID column as rownames
rownames(pub_rosmap_mgp_estimations) <- pub_rosmap_mgp_estimations$specimenID
pub_rosmap_mgp_estimations$specimenID <- NULL
rownames(sn_proportions_forpub) <- sn_proportions_forpub$specimenID
sn_proportions_forpub$specimenID <- NULL

sn_proportions_forpub$IT <- apply(sn_proportions_forpub[, c("L5.IT", "L6.IT", "L2.3.IT")], 1, sum)
sn_proportions_forpub <- sn_proportions_forpub[, !(names(sn_proportions_forpub) %in% c("L5.IT", "L6.IT", "L2.3.IT"))]

# Match up column names for sn_proportions_forpub & rosmap_mgp_estimations
setdiff(colnames(sn_proportions_forpub), colnames(pub_rosmap_mgp_estimations))
setdiff(colnames(pub_rosmap_mgp_estimations), colnames(sn_proportions_forpub))
common_subclass <- intersect(colnames(sn_proportions_forpub), colnames(pub_rosmap_mgp_estimations))
common_specimen <- intersect(rownames(sn_proportions_forpub), rownames(pub_rosmap_mgp_estimations))


pub_rosmap_mgp_estimations_common <- pub_rosmap_mgp_estimations[common_specimen, common_subclass]
sn_proportions_forpub_common <- sn_proportions_forpub[common_specimen, common_subclass]

# Using markers from Micaela's paper
pub_sn_cor <- diag(cor(pub_rosmap_mgp_estimations_common, sn_proportions_forpub_common))
pub_sn_cor <- pub_sn_cor[!is.na(pub_sn_cor)]



rosmap_mgp_estimations <- readRDS("/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/rosmap_mgp_estimations.rds") %>%
  dplyr::rename(Micro = Micro.PVM)

# Set specimenID column as rownames
rownames(rosmap_mgp_estimations) <- rosmap_mgp_estimations$specimenID
rosmap_mgp_estimations$specimenID <- NULL
rownames(sn_proportions) <- sn_proportions$specimenID
sn_proportions$specimenID <- NULL

# Match up column names for sn_proportions & rosmap_mgp_estimations
setdiff(colnames(sn_proportions), colnames(rosmap_mgp_estimations))
setdiff(colnames(rosmap_mgp_estimations), colnames(sn_proportions))
common_subclass <- intersect(colnames(sn_proportions), colnames(rosmap_mgp_estimations))
common_specimen <- intersect(rownames(sn_proportions), rownames(rosmap_mgp_estimations))

rosmap_mgp_estimations_common <- rosmap_mgp_estimations[common_specimen, common_subclass]
sn_proportions_common <- sn_proportions[common_specimen, common_subclass]

# Using markers I generated
mgp_sn_cor <- diag(cor(rosmap_mgp_estimations_common, sn_proportions_common))
mgp_sn_cor <- mgp_sn_cor[!is.na(mgp_sn_cor)]




# dtangle
rosmap_dt_estimations <- readRDS("/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/rosmap_dt_estimations.rds") %>%
  as.data.frame() %>%
  dplyr::rename(Micro = Micro.PVM) 

# Match up column names for sn_proportions & rosmap_mgp_estimations
setdiff(colnames(sn_proportions), colnames(rosmap_dt_estimations))
setdiff(colnames(rosmap_dt_estimations), colnames(sn_proportions))
common_subclass <- intersect(colnames(sn_proportions), colnames(rosmap_dt_estimations))
common_specimen <- intersect(rownames(sn_proportions), rownames(rosmap_dt_estimations))

rosmap_dt_estimations_common <- rosmap_dt_estimations[common_specimen, common_subclass]
sn_proportions_common <- sn_proportions[common_specimen, common_subclass]

# Using markers I generated
dt_sn_cor <- diag(cor(rosmap_dt_estimations_common, sn_proportions_common))
dt_sn_cor <- dt_sn_cor[!is.na(dt_sn_cor)]



# Harmonizing subclass names across the three methods
setdiff(names(pub_sn_cor), names(mgp_sn_cor))
setdiff(names(mgp_sn_cor), names(pub_sn_cor))
setdiff(names(dt_sn_cor), names(mgp_sn_cor))
setdiff(names(dt_sn_cor), names(pub_sn_cor))



all_names <- union(names(pub_sn_cor), union(names(mgp_sn_cor), names(dt_sn_cor)))

pub_sn_cor_names <- all_names[!all_names %in% names(pub_sn_cor)]
# Create vector of NA values with names
na_values <- rep(NA, length(pub_sn_cor_names))
names(na_values) <- pub_sn_cor_names
# Concatenate new NA values to existing vector
pub_sn_cor <- c(pub_sn_cor, na_values)

mgp_sn_cor_names <- all_names[!all_names %in% names(mgp_sn_cor)]
# Create vector of NA values with names
na_values <- rep(NA, length(mgp_sn_cor_names))
names(na_values) <- mgp_sn_cor_names
# Concatenate new NA values to existing vector
mgp_sn_cor <- c(mgp_sn_cor, na_values)

dt_sn_cor_names <- all_names[!all_names %in% names(dt_sn_cor)]
# Create vector of NA values with names
na_values <- rep(NA, length(dt_sn_cor_names))
names(na_values) <- dt_sn_cor_names
# Concatenate new NA values to existing vector
dt_sn_cor <- c(dt_sn_cor, na_values)



# Combine the vectors into a data frame
comb_df <- data.frame(mgp_sn_cor, pub_sn_cor, dt_sn_cor)
comb_df <- comb_df[!rownames(comb_df) %in% c("L5.ET", "L4.IT", "OPC", "IT", "Sst.Chodl"),]
comb_df <- rownames_to_column(comb_df, var = "cell_type")

comb_df_long <- pivot_longer(comb_df, cols = -cell_type, names_to = "method", values_to = "corr")

ggplot(comb_df_long, aes(x = cell_type, y = corr, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(. ~ method, scales = "free_x", space = "free_x") +
  labs(x = "Cell type", y = "Pearson's correlation") +
  theme_minimal()

# create the barplot with grouped bars
ggplot(comb_df_long, aes(x = cell_type, y = corr, fill = method, group = cell_type)) +
  geom_bar(stat = "identity", position="dodge2") +
  labs(x = "Cell Type", y = "Pearson's Correlation") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"))  # specify colors for each method
