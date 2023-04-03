
analysis_dir <- "/Users/fu/Library/CloudStorage/GoogleDrive-fu@broadinstitute.org/Shared drives/GPP Cloud /R&D/People/Dany/RNAseq Analysis/TP53 Base Editing/"
reads_file_dir <- "/Users/fu/Library/CloudStorage/GoogleDrive-fu@broadinstitute.org/Shared drives/GPP RNA Seq/2023/TP53/Terra/RSEM_Results/All-reads/"
annot <- read.csv("../TP53-BE-tiling-annotation.csv")

read_files <- list.files(path=reads_file_dir, full.names=TRUE)
reads_list <- lapply(read_files, function(x) {
  r <- fread(x, select = c("gene_id", "expected_count"))
  file_name <- unlist(strsplit(x, "\\."))[2]
  # Change the expected_count column to the sample name
  setnames(r, c("gene_id", sub(".*TP53-A549-", "", file_name)))
})

reads_df <- Reduce(function(x, y) merge(x, y, by = "gene_id", all=TRUE), reads_list)
reads_df<- as.data.frame(reads_df)
rownames(reads_df) <- reads_df[,1]
reads_df[,1] <- NULL

# Make DGE object
reads.dge <- DGEList(counts=reads_df, samples = annot)

# Remove lowly expressed genes
reads.dge$samples$group <- as.factor(c(annot$Sample))
keep.exprs.reads <- filterByExpr(reads.dge, group = reads.dge$samples$group)
reads.dge <- reads.dge[keep.exprs.reads,, keep.lib.sizes=FALSE]

# Normalize by TMM
reads.dge <- calcNormFactors(reads.dge, method = "TMM")
reads.dge.lcpm <- cpm(reads.dge, log=TRUE)

############################
# MULTIDIMENSIONAL SCALING #
############################
# Set levels to colors
sample.groups <- reads.dge$samples$group
levels(sample.groups) <-  brewer.pal(nlevels(sample.groups), "Dark2")
sample.groups <- as.character(sample.groups)

png(file = "MDS-by-replicate.png")
plotMDS(reads.dge.lcpm, col=sample.groups, cex=0.7, xlim=c(-1.5, 1.5))
dev.off()

# Plot by Exon targeted 
exon.groups <- as.factor(c(annot$Target.Exon))
levels(exon.groups) <-  brewer.pal(nlevels(exon.groups), "Dark2")
exon.groups <- as.character(exon.groups)
plotMDS(reads.keep.norm, col=exon.groups)

# Plot by treatment
cond.groups <- as.factor(c(annot$Condition))
levels(cond.groups) <-  brewer.pal(nlevels(cond.groups), "Dark2")
cond.groups <- as.character(cond.groups)
plotMDS(reads.keep.norm, col=cond.groups)

# 569-Nutlin-Reps looks weird, so do a correlation analysis

#######################
# PEARSON CORRELATION #
#######################
png(file = "569-pearson-heatmap.png")
all_569_counts <- reads.keep.norm[,c("RDA569-untreated-RepA", "RDA569-untreated-RepB", "RDA569-untreated-RepC",
                                            "RDA569-Nutlin-RepA", "RDA569-Nutlin-RepB", "RDA569-Nutlin-RepC")]
cor_569_all <- cor(all_569_counts)
pheatmap(cor_569_all)
dev.off()

png(file = "692-pearson-heatmap.png")
all_692_counts <- reads.keep.norm[,c("RDA692-untreated-RepA", "RDA692-untreated-RepB", "RDA692-untreated-RepC",
                                            "RDA692-Nutlin-RepA", "RDA692-Nutlin-RepB", "RDA692-Nutlin-RepC")]
cor_692_all <- cor(all_692_counts)
pheatmap(cor_692_all)
dev.off()

png(file = "429-pearson-heatmap.png")
all_429_counts <- reads.keep.norm[,c("RDA429-untreated-RepA", "RDA429-untreated-RepB", "RDA429-untreated-RepC",
                                            "RDA429-Nutlin-RepA", "RDA429-Nutlin-RepB", "RDA429-Nutlin-RepC")]
cor_429_all <- cor(all_429_counts)
pheatmap(cor_429_all)
dev.off()

png(file = "pearson-heatmap-all.png")
cor_all <- cor(reads.keep.norm)
pheatmap(cor_all)
dev.off()






##########################
# TP53 & MDM2 EXPRESSION #
##########################
tp53.id <- "ENSG00000141510" # ENSG00000141510.17
mdm2.id <- "ENSG00000135679"

reads.tp53 <- reads.dge.lcpm[grepl(tp53.id, rownames(reads.dge.lcpm)),]
reads.tp53.long <- data.frame("reads" = reads.tp53)
reads.tp53.long$transcript_id <- rownames(reads.tp53.long)
reads.tp53.long$group <- reads.dge$samples$group
reads.tp53.long$condition <- annot$Condition
reads.tp53.long <- melt(setDT(reads.tp53.long), 
                        id.vars = c("transcript_id", "group", "condition"), 
                        variable.name = "reads")

reads.mdm2 <- reads.dge.lcpm[grepl(mdm2.id, rownames(reads.dge.lcpm)),]
reads.mdm2.long <- data.frame("reads" = reads.mdm2)
reads.mdm2.long$transcript_id <- rownames(reads.mdm2.long)
reads.mdm2.long$group <- reads.dge$samples$group
reads.mdm2.long$condition <- annot$Condition
reads.mdm2.long <- melt(setDT(reads.mdm2.long), 
                        id.vars = c("transcript_id", "group", "condition"), 
                        variable.name = "reads")


pdf(file = paste(analysis_dir, "boxplot.pdf"))
tp53.bp <- boxplot(reads.tp53.long, title="TP53 Expression")
mdm2.bp <- boxplot(reads.mdm2.long, title="MDM2 Expression")
ggarrange(tp53.bp, mdm2.bp, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off() 



# Remove RDA-569Nutlin-RepA
rda569_nutlinA.rownum <- which(rownames(reads.tp53.long) == "RDA569-Nutlin-RepA") 
reads.tp53.long.filtered <- reads.tp53.long[-rda569_nutlinA.rownum,]
png(file = paste(analysis_dir, "tp53-boxplot-outlier-removed.png"))
boxplot(reads.tp53.long.filtered)
dev.off()

boxplot = function(data, title) {
  return (
    ggplot(data, aes(x=group, y=value, palette = "BuGn",fill=condition)) + 
      geom_boxplot() + 
      ylab("log(cpm)") + xlab ("Samples") + labs(title = title) + 
      theme(legend.title = element_text(size = 15), 
              legend.text = element_text(size = 15)) + 
      scale_y_continuous(limits = c(4.9, 9.9)) + 
      scale_x_discrete(guide = guide_axis(angle = 90),
                       labels= c("Control", "Control", "SG1", "SG1","SG9","SG9"), 
                       limits=c("untreated_ctrl","treated_ctrl",
                                "untreated_ex1", "treated_ex1",
                                "untreated_ex9", "treated_ex9")) + 
      stat_compare_means(label = "p.signif",method = "t.test", 
                         comparisons = list(c("treated_ctrl", "untreated_ctrl"),
                                            c("treated_ex1", "untreated_ex1"),
                                            c("treated_ex9", "untreated_ex9")))
  )
}


