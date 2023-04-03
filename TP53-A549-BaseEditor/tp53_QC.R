library(limma)
library(data.table)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(edgeR)

qc_file_dir <- "/Users/fu/Library/CloudStorage/GoogleDrive-fu@broadinstitute.org/Shared drives/GPP RNA Seq/2023/TP53/Terra/RNASeQC_Results/"
qc_files <- list.files(path=qc_file_dir, full.names=TRUE)
# read each file using the first column as the header, then transpose
qc_df_list <- lapply(qc_files, function(x) {
  t(read.csv(x, row.names=1, sep = '\t'))
})
# Merge all the dataframes together
qc_df <- as.data.frame(do.call(rbind, qc_df_list))

#####################
# Outlier histogram #
#####################
# Remove columns where all the values are the same
qc_df.unique <- qc_df[vapply(qc_df, function(x) length(unique(x)) > 1, logical(1L))]
qc_tally <- data.frame(matrix(NA, nrow = nrow(qc_df.unique), ncol = ncol(qc_df.unique)))
rownames(qc_tally) <- rownames(qc_df.unique)
colnames(qc_tally) <- colnames(qc_df.unique)
for(i in 1:ncol(qc_df.unique)) { 
  # Find mean and SD for every column (stat)
  sd.qc <- sd(qc_df.unique[,i])
  mean.qc <- mean(qc_df.unique[,i])
  
  # For each row, find if sample exceed 2sd above or below the mean
  qc_col <- qc_df.unique[, i]
  qc_tally[,i] <- (qc_col > (mean.qc + 2 * sd.qc) | qc_col < (mean.qc - 2 * sd.qc))
}

# simple bargraph 
barplot(rowSums(qc_tally),las=2)

# histogram
qc.outlier_sums <- data.frame(count=rowSums(qc_tally))
qc.outlier_histo <- ggplot(qc.outlier_sums, aes(x=count)) + 
  geom_histogram() +
  xlab("Tally of outliers")
ggsave("outlier_histo.pdf", qc.outlier_histo, width=10)

##################################
# Violin plot of number of genes #
##################################
df.genes_detected <- qc_df["Genes Detected"]
df.genes_detected$Samples <- rownames(qc_df)
df.mapped_rates.long <- melt(setDT(df.genes_detected), id.vars = "Samples", measure.vars = "Genes Detected")
numGenes.vplot <- ggplot(df.mapped_rates.long, aes(x=variable, y=value)) + geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 2), geom="pointrange", color="red") +
  labs(x = "Samples", y = "Number of Genes Detected") + 
  theme(axis.text.x = element_blank())
# numGenes.vplot + stat_summary(fun.y=median, geom="point", size=2, color="red")
ggsave("numGenes_vplot.pdf", numGenes.vplot)

####################################
# Violin plot of High Quality Rate #
####################################
df.high_quality_rate <- qc_df["High Quality Rate"]
df.high_quality_rate$Samples <- rownames(qc_df)
df.high_quality_rate.long <- melt(setDT(df.high_quality_rate), id.vars = "Samples", measure.vars = "High Quality Rate")
highQualityRate.vplot <- ggplot(df.high_quality_rate.long, aes(x=variable, y=value)) + geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 2), geom="pointrange", color="red") +
  labs(x = "Samples", y = "High Quality Rate") + 
  theme(axis.text.x = element_blank())
# numGenes.vplot + stat_summary(fun.y=median, geom="point", size=2, color="red")
ggsave("highQualityGenes_vplot.pdf", highQualityRate.vplot)

sd(df.high_quality_rate$`High Quality Rate`)
mean(df.high_quality_rate$`High Quality Rate`)

################
# RATE METRICS #
################
# Heatmap
all_rate_cols <- c("Mapping Rate", "Exonic Rate", "Intronic Rate", "Intergenic Rate", 
                  "Ambiguous Alignment Rate", "rRNA Rate")
heatmap(qc_df[,all_rate_cols])

# Stacked bar graph
mappedRateCols <- c("Ambiguous Alignment Rate","Intergenic Rate","Intronic Rate","Exonic Rate")
df.mapped_rates <- qc_df[,mappedRateCols]
df.mapped_rates$Samples <- rownames(df.mapped_rates)
# Convert wide to long for stacked bar plot
df.mapped_rates.long <- melt(setDT(df.mapped_rates), id.vars = "Samples", measure.vars = mappedRateCols)
alignmentRate.stacked <- ggplot(df.mapped_rates.long,aes(fill=variable, x=value, y=Samples )) + 
  geom_bar(position="stack", stat="identity")+theme_classic()
ggsave("QC/alignmentRate_stacked.pdf", alignmentRate.stacked, width=10)



