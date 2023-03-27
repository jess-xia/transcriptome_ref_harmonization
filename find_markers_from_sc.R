library(dplyr)
library(Seurat)
library(patchwork)
library(MAST)
library(org.Hs.eg.db)
library(dplyr)

# Load the single cell meta data and counts
sc_meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/Isoform_project/Interspecies_comparison/human/metadata.csv")
sc_exon <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/Isoform_project/Interspecies_comparison/human/human_bams/M1/exon_counts.csv")
sc_intron <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/Isoform_project/Interspecies_comparison/human/human_bams/M1/intron_counts.csv")

# Metadata
sc_meta <- sc_meta[,-1] # Remove first column X in dataframe

table(sc_meta$outlier_call, sc_meta$region_label) # Check for outliers
new_metadata_filtered <- new_metadata[new_metadata$outlier_call == "False",] # Keep only non-outliers

table(new_metadata_filtered$subclass_label) #Look at sample counts for each subclass
table(new_metadata_filtered[,c("subclass_label", "region_label")], useNA = "ifany")
row.names(new_metadata_filtered) <- new_metadata_filtered$sample_name #setting row names to sample names

# Combine the exon and intron counts into one dataframe
sc_sum <- rbind(sc_exon, sc_intron)
sc_sum <- sc_sum[,-1] # Remove first column X in dataframe
# sc_allcounts_raw <- aggregate(x = sc_sum[ , colnames(sc_sum) != "geneid"], # Sum by group
#           by = list(sc_sum$geneid),
#           FUN = sum)
# saveRDS(sc_allcounts, "~/transcriptome-ref-harmonization/sc_allcounts_raw.rds")
sc_allcounts_raw <- readRDS("~/transcriptome-ref-harmonization/sc_allcounts_raw.rds")
sc_allcounts <- column_to_rownames(sc_allcounts_raw, "Group.1") # Set first column as rownames, needed to convert to matrix
sc_allcounts <- sc_allcounts[,rownames(new_metadata_filtered)] #filter for the cells we specified previously (removing outliers)


# Convert dataframe to matrix 
new_count_matrix <- as.matrix(sc_allcounts, sparse = TRUE) 


# Create seurat object from counts and metadata
humanM1_seuobj <- CreateSeuratObject(counts = new_count_matrix, meta.data = new_metadata_filtered) 
humanM1_seuobj <- NormalizeData(humanM1_seuobj, normalization.method = "LogNormalize", scale.factor = 1000000)

Idents(humanM1_seuobj) <- "subclass_label" # set active identity
head(Idents(humanM1_seuobj)) # see our active identities of cells

all_group_markers <- FindAllMarkers(object = humanM1_seuobj, 
                                    logfc.threshold = 2.5, # only test a gene if there is a log-fold-change difference above this threshold
                                    min.pct = .35, # only test a gene if the % of cells in the group of interest that express it is > this threshold
                                    only.pos = T, # only return positive markers (enriched in the cell group of interest)
                                    test.use = "roc") # the statistical test to use
saveRDS(all_group_markers, "~/transcriptome-ref-harmonization/M1_all_group_markers.rds")



sonny_markers <- read_csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'))
marker_data <- readRDS("/external/rprshnas01/kcni/jxia/transcriptome-ref-harmonization/M1_all_group_markers.rds")

annots <- select(org.Hs.eg.db, keys=marker_data$gene, 
                 columns="SYMBOL", keytype="ENSEMBL")
# Try again after dplyr is fixed
result <- left_join(marker_data, annots, by = c("gene"="ENSEMBL"))
