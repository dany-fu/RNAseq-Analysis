library(biomaRt)

isoforms_file_dir <- "/Users/fu/Library/CloudStorage/GoogleDrive-fu@broadinstitute.org/Shared drives/GPP RNA Seq/2023/TP53/Terra/RSEM_Results/Isoforms/"
isoform_files <- list.files(path=isoforms_file_dir, full.names=TRUE)

# Merge the transcript reads into one dataframe
isoforms_list <- lapply(isoform_files, function(x) {
  r <- fread(x, select = c("transcript_id", "expected_count"))
  file_name <- unlist(strsplit(x, "\\."))[2]
  # Replace the expected_count column name to the sample name
  setnames(r, c("transcript_id", sub(".*TP53-A549-", "", file_name)))
})
isoforms_merged <- Reduce(function(x, y) merge(x, y, by = "transcript_id", all=TRUE), isoforms_list)
isoforms_df<- as.data.frame(isoforms_merged)
# Turn the first column 'transcript_id' to rownames
rownames(isoforms_df) <- isoforms_df[,1]
isoforms_df[,1] <- NULL
# sort the column names by alphabetical order
isoforms_df <- isoforms_df[ , sort(names(isoforms_df))]

# Make DGE object
isoforms.dge <- DGEList(counts=isoforms_df)

# Remove lowly expressed genes
isoforms.dge$samples$group <- as.factor(c(annot$Sample))
keep.exprs.isoform <- filterByExpr(isoforms.dge, group = isoforms.dge$samples$group)
isoforms.dge.keep <- isoforms.dge[keep.exprs.isoform,, keep.lib.sizes=FALSE]

# Normalize by TMM
isoforms.dge.keep <- calcNormFactors(isoforms.dge.keep, method = "TMM")
isoforms_dge.keep.tmm <- cpm(isoforms.dge.keep, log=F)

# Get all TP53 transcript IDs from Ensembl
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
tp53.info <- getBM(attributes=c('ensembl_transcript_id_version','ensembl_gene_id', 'ensembl_transcript_id'),filters = 'ensembl_gene_id', values = tp53.id, mart = ensembl)

# Remove the version number from the transcript ID
# isoforms_df$transcript_id <- sapply(isoforms_df$transcript_id, function(x) {
#   unlist(strsplit(x, "\\."))[1]
# })

# Subset to only the TP53 isoforms
match_pattern <- paste(tp53.info$ensembl_transcript_id,collapse="|")
tp53.isoforms <- isoforms_dge.keep.tmm[grepl(match_pattern, rownames(isoforms_dge.keep.tmm)),]
pheatmap(tp53.isoforms, clustering_method = "ward.D2")

# Look at TP53 isoform raw reads
tp53.isoforms.raw <- isoforms_df[grepl(match_pattern, rownames(isoforms_df)),]

# Stacked bargraph of the isoforms
tp53.isoforms.long <- data.frame(t(tp53.isoforms))
tp53.isoforms.long$transcript_id <- rownames(tp53.isoforms.long)
tp53.isoforms.long <- melt(setDT(tp53.isoforms.long), id.vars="transcript_id", measure.vars=rownames(tp53.isoforms))

untreated429.repA <- tp53.isoforms[, 6]
untreated429.repA <- untreated429.repA[order(untreated429.repA)]
tp53.level_order <- names(untreated429.repA)

png(file = paste(analysis_dir, "isoform-stacked-bargraph.png"), width = 600)
ggplot(tp53.isoforms.long,aes(x=value, y=transcript_id,
                              fill = factor(variable, levels=tp53.level_order))) + 
  geom_bar(position="fill", stat="identity")+theme_classic() # fill for percentage
dev.off()

