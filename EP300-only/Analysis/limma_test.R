# https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#differential-expression-analysis
library(limma)
library(Glimma)
library(edgeR)
library(pheatmap)
library(gplots)
library(viridis)
library(SummarizedExperiment)
library(biomaRt)
library(RColorBrewer)
library(ggpubr)

setwd("/Users/danyfu/Dropbox/CompBioMed/RNA_Seq/JnJInvitro/EP300-only/Analysis/")
ep300.exp<-readRDS(file="JnJ_invitro_ep300__Gene_Expression.rds")
ep300.raw_counts <- assay(ep300.exp)
rna.all.annot<-read.delim(file="annotation.txt",sep=",")
# dge.all <- DGEList(counts=ep300.raw_counts, samples=rna.all.annot)
# dge.all.keep <- removeGenes(dge.all)
# which(rownames(dge.all.keep) == "ENSG00000100393") 
# dge.all.keep <- calcNormFactors(dge.all.keep, method = "TMM")
# dge.all.keep.tmm <- cpm(dge.all.keep, log=TRUE) 
# makeMDSPlot(dge.all.keep.tmm, dge.all.keep)



###############################################
# 5.2 Removing genes that are lowly expressed #
###############################################
# Determine which genes have sufficiently large counts to be retained, group by cell line
removeGenes = function (dge, grouping) {
  keep.exprs <- filterByExpr(dge, group = grouping)
  # false means library size will be recalculated to be the sum of the counts left in the rows 
  return(dge[keep.exprs,, keep.lib.sizes=FALSE])
}

#################################################
# 5.3 Normalising gene expression distributions #
#################################################
# normalisation by the method of trimmed mean of M-values (TMM), scaling factor for the library sizes
# pass in the filtered DEG

##########################################
# 5.4 Unsupervised clustering of samples #
##########################################
# multi-dimensional scaling (MDS) plot
# pass in DGE after filtering and normalization
makeMDSPlot = function (lcpm, dge){
  par(mfrow=c(1,2))
  col.group <- as.factor(c(dge$samples$Sample))
  levels(col.group) <- brewer.pal(nlevels(col.group), "Paired") #Paired
  col.group <- as.character(col.group)
  colors = c("#B3B3B3", "#B3B3B3","#B3B3B3",
             "#A6D854", "#A6D854", "#A6D854",
             "#FFD92F", "#FFD92F", "#FFD92F",
             "#E5C494", "#E5C494", "#E5C494",
             "#3399FF", "#3399FF","#3399FF",
             "#FC8D62", "#FC8D62","#FC8D62",
             "#CC0099", "#CC0099","#CC0099",
             "#000033", "#000033", "#000033")
  return(plotMDS(lcpm, labels=dge$samples$SampleID, col=colors,
                 cex=0.7, xlim=c(-4, 4)))
}



######################################
# 6 Differential expression analysis #
######################################
# CellLine <- dge.interest.withAnnot.keep.normalized$samples$CellLine # add this if looking at dge.all.withAnnot

# No intercept 
#design_0 <- model.matrix(~0+Perturbation+CellLine)
#contr.matrix <- makeContrasts(ControlvsEP300 = PerturbationControl-PerturbationEP300, levels = colnames(design))

transformAndFit = function(dge) {
  Perturbation <- dge$samples$Perturbation
  #  intercept is the first parameter
  design_1 <- model.matrix(~Perturbation)
  v<- voom(dge, design_1, plot=FALSE)
  dupCor <- duplicateCorrelation(v,design_1,block=dge$samples$Sample)
  fitRan <- lmFit(v,design_1,block=dge$samples$Sample,correlation=dupCor$consensus)
  return(eBayes(fitRan))
}

###########################################################################
# 6.6 Useful graphical representations of differential expression results #
###########################################################################

# convert ensemble IDs to HGNC gene names
convertEnsIdToGeneName = function(ensembleId.list){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genelist <- getBM(filters=c("ensembl_gene_id"), attributes= c("ensembl_gene_id","hgnc_symbol"),values=ensembleId.list,mart=mart)
  return(genelist[match(ensembleId.list, genelist$ensembl_gene_id),])
}

# # heatmap of the top genes
# colorPalette <- colorpanel(1000,"blue","white","red")
# col.sample <- as.factor(c(dge.interest.withAnnot.keep.normalized$samples$Sample))
# heatmap.2(lcpm[i,], scale="row",cexCol=0.75,labCol=col.sample, labRow=FALSE, 
#           col=colorPalette, trace="none", density.info="none", 
#           margin=c(8,6), lhei=c(2,10), dendrogram="column",
#           distfun = function(x) dist(x, method="euclidean"),
#           hclustfun = function(x) hclust(x, method="ward.D2"))
# 
# # heatmap of the top 50 down regulated
# ep300Genes.down <- subset(ep300Genes, logFC <0)
# top50.down <- ep300Genes.down[order(ep300Genes.down$adj.P.Val),][1:50,]
# top50.down.geneNames <- rownames(top50.down)
# d <- which(rownames.voom.all %in% top50.down.geneNames)
# heatmap.2(lcpm[d,], scale="row",cexRow=0.5,cexCol=0.75,
#           labRow=ens2$hgnc_symbol[d], labCol=col.group, 
#           col=colorPalette, trace="none", density.info="none", 
#           margin=c(8,6), lhei=c(2,10), dendrogram="column")



