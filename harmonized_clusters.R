library(dplyr)
library(readr)
library(vctrs)
library(tidyr)

# Import sn proportions data, convert to long, and match up the microglial name
sn_proportions_raw <- read_csv("/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/rosmap_single_nuc_proportions.csv") %>%
  subset(select = -X1)
sn_proportions_long <- gather(sn_proportions_raw, subclass, estimation, Astro:Oligo, factor_key=FALSE) %>%
  filter(!subclass %in% c("Lamp5 Lhx6", "Chandelier", "Pax6", "L4 IT" )) %>% # Harmonize subclass names 
  rename(sn_subclass = subclass) %>%
  subset()
sn_proportions_long$subclass[sn_proportions_long$subclass == "Micro"] <- "Micro-PVM" # Match up microglia subclass name 

# Load Azimuth cell type proportion predictions
azimuth <- read_tsv("/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/azimuth_pred.tsv") %>%
  filter(predicted.subclass != "VLMC") %>% # Harmonize subclass names
  rename(azimuth_subclass = predicted.subclass)
# Compare difference between the two list of subclasses
sn_types <- unique(sn_proportions_long$subclass)
azimuth_types <- unique(azimuth$predicted.subclass)
setdiff(sn_types, azimuth_types)
setdiff(azimuth_types, sn_types)

setdiff(sn_types, colnames(rosmap_mgp_estimations))

combined_sst <- merge(sn_sst, mgp_sst, by = "specimenID")


colnames(rosmap_mgp_estimations)
unique(sn_proportions_long$sn_subclass)
