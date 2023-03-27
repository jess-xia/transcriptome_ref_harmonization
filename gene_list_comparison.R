# Import bulk and sc RNA-seq data 
# Note: Number after the decimal point in ensg id is the verion, can be ignored
bulk <- read.delim("/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/raw_gene_counts/gene_quantification/rosmap_batch1_stranded/ROSMAP_batch1_gene_all_counts_matrix.txt")
sc_meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/Isoform_project/Interspecies_comparison/human/metadata.csv")
# sc_intron & sc_exon have the same number of observations, same gene list
sc_exon <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/Isoform_project/Interspecies_comparison/human/human_bams/M1/exon_counts.csv")
sc_intron <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/Isoform_project/Interspecies_comparison/human/human_bams/M1/intron_counts.csv")

 
# Make new column with just ENSG ID, no version
bulk$gene_id <- sub("\\..*", "", bulk$feature)
# Check if single cell data has any versions specified
sum(sapply(sc_exon$geneid, grepl, pattern = "//."))

# Summary for shared vs unique between the two datasets
total_count <- length(union(bulk$gene_id, sc_exon$geneid))
percent_shared <- length(intersect(bulk$gene_id, sc_exon$geneid))/total_count * 100
percent_bulk_unique <- length(setdiff(bulk$gene_id, sc_exon$geneid))/total_count * 100
percent_sc_unique <- length(setdiff(sc_exon$geneid, bulk$gene_id))/total_count * 100

