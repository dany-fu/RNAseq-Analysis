library(GSVA)
library(data.table)
library(nlme)
library(pheatmap)
library(RColorBrewer)
library(SummarizedExperiment)
library(ggpubr)

setwd("/rprojectnb2/pcga/JnJ_invitro/Analysis/GSVA")
# Using DEGs from NL20 c11 and c18
nl20.degs<-readRDS(file="nl20-best-degs.rds") #1191 genes
ep300.ensId <-"ENSG00000100393"

ep300Genes.down <- rownames(subset(nl20.degs, logFC <0)) #847
ep300Genes.up <- rownames(subset(nl20.degs, logFC >0)) #344

########
# CCLE #
########
ccle.info<-read.delim(file="sample_info.csv",sep=",")
# ccle.lung <-ccle.info[ccle.info$primary_disease=="Lung Cancer",]
# ccle.lung.lusc <- ccle.lung[ccle.lung$Subtype=="Non-Small Cell Lung Cancer (NSCLC), Squamous Cell Carcinoma",]
ccle.info.lusc <- ccle.info[ccle.info$lineage_sub_subtype=="NSCLC_squamous",]

ccle.exp<-read.delim(file="CCLE_RNAseq_reads.csv",sep=",")
ccle.lusc <- ccle.exp[ccle.exp$X %in% ccle.info.lusc$DepMap_ID,]
# only keep the ENS ids in the gene names
# colnames(ccle.lusc) <- gsub("\\..*","", colnames(ccle.lusc)) # for CCLE_expression_full.csv
colnames(ccle.lusc) <- gsub(".*\\.(.+)\\.", "\\1", colnames(ccle.lusc))
# transpose the matrix so that rows are gene names
ccle.lusc <- setNames(data.frame(t(ccle.lusc[,-1])), ccle.lusc[,1])

# write just the LUSC lines for easier access later
saveRDS(data.matrix(ccle.lusc), file="CCLE_LUSC_raw_reads.rds")
ccle.lusc<-readRDS(file="CCLE_LUSC_raw_reads.rds")

# Do all the normalization and log TMM stuff
ccle.lusc.dge <- DGEList(counts=ccle.lusc)
ccle.lusc.keep <- removeGenes(ccle.lusc.dge, colnames(ccle.lusc.dge))
ccle.lusc.keep <- calcNormFactors(ccle.lusc.keep, method = "TMM")
ccle.lusc.keep.tmm <- cpm(ccle.lusc.keep, log=TRUE) 

# write the normalized LUSC lines
saveRDS(ccle.lusc.keep.tmm,file="CCLE_norm_reads.rds")

#Get the common set of all genes between the EP300 data and CCLE
ccle.lusc.keep.tmm<-readRDS(file="CCLE_norm_reads.rds")
# get gene names for nl20 DEGs
#degNames <- convertEnsIdToGeneName(rownames(nl20.degs))
#degNames <- degNames[!degNames$hgnc_symbol == "", ] # remove genes that do not have a HGNC name
common.table.ccle<-nl20.degs[rownames(nl20.degs) %in% rownames(ccle.lusc.keep.tmm),] #1191

#Create a list of DE invitro genes.  One list of sign. up-reg. genes, one of sign. dn-reg, and one of all genes.
common.list.ccle<-c(list(rownames(common.table.ccle)[which(common.table.ccle$t>0)]),
               list(rownames(common.table.ccle)[which(common.table.ccle$t<0)]),
               list(rownames(common.table.ccle)))
names(common.list.ccle)<-c("UP","DN","ALL")

#Run GSVA on CCLE LUSC data
gsva.ccle<-gsva(data.matrix(ccle.lusc.keep.tmm), common.list.ccle) # 27 cell lines
saveRDS(gsva.ccle,file="GSVA_results_nl20_ccle.rds")

# Plot GSVA results on heatmap
gsva.ccle<-readRDS(file="GSVA_results_nl20_ccle.rds")
pheatmap(gsva.ccle, clustering_method = "ward.D2")

#Append EP300 expression, transpose, and create dataframe for the scatter plot
gsva.ccle.ep300 <- data.frame(t(rbind(gsva.ccle, EP300=ccle.lusc.keep.tmm[ep300.ensId,])))
saveRDS(gsva.ccle.ep300,file="GSVA_results_nl20_ccle_with_EP300.rds")

gsva.ccle.ep300<-readRDS(file="GSVA_results_nl20_ccle_with_EP300.rds")
# Plot scatter graph with Pearson correlation between EP300 expression and GSVA scores
pdf(file="GSVA_plots_CCLE_NL20.pdf")
fillColors <- c("bisque","#b3ecec","#e7d1ff")
par(mfrow=c(1,3))
for(i in 1:3){
  print(ggscatter(gsva.ccle.ep300, x = "EP300", y = colnames(gsva.ccle.ep300)[i], 
                  title=paste("CCLE LUSC EP300",colnames(gsva.ccle.ep300)[i], sep=" "),
                  add = "reg.line", conf.int = TRUE, add.params = list(fill=fillColors[i]),
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "EP300 Expression (log CPM)", ylab = "GSVA Score"))
}
dev.off()

########
# PCGA #
########
pcga<-readRDS(file="/restricted/projectnb/pcga/ANALYSIS_NEXTFLOW/edgeRobj/cpm_bx_all.rds")
common.table.pcga <- nl20.degs[rownames(nl20.degs) %in% rownames(pcga), ] # 990 genes
common.list.pcga<-c(list(rownames(common.table.pcga)[which(common.table.pcga$t>0)]),
               list(rownames(common.table.pcga)[which(common.table.pcga$t<0)]),
               list(rownames(common.table.pcga)))
names(common.list.pcga)<-c("UP","DN","ALL")

#Run GSVA on PCGA data
gsva.pcga<-gsva(pcga, common.list.pcga)
saveRDS(gsva.pcga,file="GSVA_results_nl20_pcga.rds")

# Plot GSVA results on heatmap
#gsva.pcga<-readRDS(file="GSVA_results_nl20_pcga.rds")
pheatmap(gsva.pcga, show_colnames=F, clustering_method = "ward.D2")

#Append EP300 expression and transpose
gsva.pcga.ep300 <- t(rbind(gsva.pcga, EP300=pcga[ep300.ensId,]))
saveRDS(gsva.pcga.ep300,file="GSVA_results_nl20_pcga_with_EP300.rds")

# Plot scatter graph with Pearson correlation between EP300 expression and GSVA scores
gsva.pcga.ep300<-readRDS(file="GSVA_results_nl20_pcga_with_EP300.rds")
fillColors <- c("bisque","#b3ecec","#e7d1ff")
pdf(file="GSVA_scatter_plots_PCGA_NL20.pdf")
plots <- list()
for(i in 1:3){
  plots[[i]] <- print(ggscatter(data.frame(gsva.pcga.ep300), x = "EP300", y = colnames(gsva.pcga.ep300)[i], 
                  title=paste("PCGA EP300",colnames(gsva.pcga.ep300)[i], sep=" "),
                  add = "reg.line", conf.int = TRUE, add.params = list(fill=fillColors[i]),
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "EP300 Expression (log CPM)", ylab = "GSVA Score"))
}
ggarrange(plotlist=plots, ncol=3, nrow=1)
dev.off()

# Get metadata for boxplots
pcga.annot<-readRDS(file="/restricted/projectnb/pcga/ANALYSIS_NEXTFLOW/subtype_associations/samples_annotation_v3.rds")
# match sample annot to what's in the counts (ie brush only)
gsva.pcga.annot <- pcga.annot[rownames(pcga.annot) %in% rownames(gsva.pcga.ep300),]
# drop all with UNK grade
gsva.pcga.annot <-subset(gsva.pcga.annot, (gsva.pcga.annot$Dysplasia_Grade != "UNK")) #291
gsva.pcga.ep300 <- gsva.pcga.ep300[rownames(gsva.pcga.ep300) %in% rownames(gsva.pcga.annot),]

grade<-as.factor(gsva.pcga.annot$Dysplasia_Grade)
grade<-relevel(grade,ref="Normal") # reorder so that normal is the first level
grade <- droplevels(grade) # drop UNK from levels
smoke<-as.factor(gsva.pcga.annot$Genomic_Smoking_Status)
tin<-as.numeric(gsva.pcga.annot$Median_TIN)
pat<-as.factor(as.character(gsva.pcga.annot$Patient))
batch<-as.factor(as.character(gsva.pcga.annot$Batch))

# Make boxplots for all 3 sets of scores & EP300 expression
y.label="GSVA Scores"
y.lim=c(-0.5, 0.5)
y.text=0.5
pcga.title="EP300 NL20 GSVA scores vs PCGA"
pdf(file="GSVA_plots_PCGA_NL20.pdf")
for(i in 1:4){
  gsva.vals<-as.numeric(gsva.pcga.ep300[,i])
  
  model.grade.1 <- lme(gsva.vals ~ grade + smoke + tin + batch, random = ~1|as.factor(pat), method="ML")
  model.grade.2 <- lme(gsva.vals ~ smoke + tin + batch, random = ~1|as.factor(pat), method="ML")
  grade.pval <- anova(model.grade.1, model.grade.2)$p[2]
  
  model.smoke.1 <- lme(gsva.vals ~ grade + smoke + tin + batch, random = ~1|as.factor(pat), method="ML")
  model.smoke.2 <- lme(gsva.vals ~ grade + tin + batch, random = ~1|as.factor(pat), method="ML")
  sm.pvals <- anova(model.smoke.1, model.smoke.2)$p[2]
  
  if(i==4){
    y.label="EP300 Expression (log CPM)"
    y.lim=c(5.5,8.5)
    y.text=8.5
  }
  graphTitle=paste(pcga.title,colnames(gsva.pcga.ep300)[i], sep=":")
  
  # Boxplot by histology
  boxplot(gsva.vals~grade,main=graphTitle,col=c("#A2CD5A", "#E0EEC7", "#FFFFFF", "#F9BDBD","#F37C7C","#EE3B3B"),
          ylab=y.label, xlab="Histology", ylim=y.lim,cex.axis=0.8)
  text(x=5, y=y.text,label=paste("Anova P value =", formatC(grade.pval,format="e",digits=2), sep=" "))
  
  # Boxplot by smoking status
  boxplot(gsva.vals~smoke, main=graphTitle, col=c("blue", "orange"),
          ylab=y.label, xlab="Smoking Status", ylim=y.lim)
  text(x=2, y=y.text,label=paste("Anova P value =", formatC(sm.pvals,format="e",digits=2), sep=" "))
}
dev.off()




#################
# Merrick et al #
#################
#This is another set of biopsy data that is published (Merrick et al./GSE114489)
merrick<-readRDS(file="/restricted/projectnb/pcga/ANALYSIS_NEXTFLOW/GSE114489/GSE114489_ExpressionSet.rds")
merrick.data<-exprs(merrick)
rownames(merrick.data)<-gsub("_at","",rownames(merrick.data))

common.table.merrick <- nl20.degs[rownames(nl20.degs) %in% rownames(merrick.data), ]
common.list.merrick<-c(list(rownames(common.table.merrick)[which(common.table.merrick$t>0)]),
                    list(rownames(common.table.merrick)[which(common.table.merrick$t<0)]),
                    list(rownames(common.table.merrick)))
names(common.list.merrick)<-c("UP","DN","ALL")

#Run GSVA on Merrick data
gsva.merrick<-gsva(merrick.data, common.list.merrick)
saveRDS(gsva.merrick,file="GSVA_results_nl20_merrick.rds")

# Plot GSVA results on heatmap
#gsva.merrick<-readRDS(file="GSVA_results_nl20_merrick.rds")
pheatmap(gsva.merrick, show_colnames=F, clustering_method = "ward.D2")

#Append EP300 expression and transpose
gsva.merrick.ep300 <- t(rbind(gsva.merrick, EP300=merrick.data[ep300.ensId,]))
saveRDS(gsva.merrick.ep300,file="GSVA_results_nl20_merrick_with_EP300.rds")

# Plot scatter graph with Pearson correlation between EP300 expression and GSVA scores
gsva.merrick.ep300<-readRDS(file="GSVA_results_nl20_merrick_with_EP300.rds")
pdf(file="GSVA_scatter_plots_Merrick_NL20.pdf")
fillColors <- c("bisque","#b3ecec","#e7d1ff")
par(mfrow=c(1,3))
for(i in 1:3){
  print(ggscatter(data.frame(gsva.merrick.ep300), x = "EP300", y = colnames(gsva.merrick.ep300)[i], 
                  title=paste("Merrick EP300",colnames(gsva.merrick.ep300)[i], sep=" "),
                  add = "reg.line", conf.int = TRUE, add.params = list(fill=fillColors[i]),
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "EP300 Expression (log CPM)", ylab = "GSVA Score"))
}
dev.off()

# Get metadata for the boxplots
gsva.merrick.ep300<-readRDS(file="GSVA_results_nl20_merrick_with_EP300.rds")
merrick.annot<-read.table(file="/restricted/projectnb/pulmseq/PCGA_Invitro/ANALYSIS/analysis101118/gse114489_annot.txt",sep="\t",header=T)
# reorder the biopsy number from the annotation to the sample name
biop.num<-c()
for(i in 1:nrow(gsva.merrick.ep300)){
  val<-strsplit(rownames(gsva.merrick.ep300)[i],"Biopsy_")
  val2<-strsplit(val[[1]][2],".CEL.gz")
  biop.num<-c(biop.num,val2[[1]][1])
}
merrick.annot<-merrick.annot[match(biop.num,merrick.annot[,1]),]

grade<-trunc(pData(merrick)$frozenDX)
grade[which(grade==7)]<-6
grade<-as.factor(grade)
smoke<-as.factor(merrick.annot$Smoke)
smoke[which(smoke=="N")]<-"F"
smoke<-as.factor(as.character(smoke))
pat<-as.factor(merrick.annot$Subject)

# Make boxplots for all 3 sets of scores & EP300 expression
y.label="GSVA Scores"
y.lim=c(-0.5, 0.5)
y.text=0.5
merrick.title="EP300 NL20 GSVA scores vs Merrick"
pdf(file="GSVA_plots_Merrick_NL20.pdf")
for(i in 1:ncol(gsva.merrick.ep300)){
  gsva.vals<-as.numeric(gsva.merrick.ep300[,i])
  
  model.grade.1 <- lme(gsva.vals ~ grade + smoke, random = ~1|as.factor(pat), method="ML")
  model.grade.2 <- lme(gsva.vals ~ smoke, random = ~1|as.factor(pat), method="ML")
  grade.pval <- anova(model.grade.1, model.grade.2)$p[2]
  print(grade.pval)
  
  model.smoke.1 <- lme(gsva.vals ~ grade + smoke, random = ~1|as.factor(pat), method="ML")
  model.smoke.2 <- lme(gsva.vals ~ grade, random = ~1|as.factor(pat), method="ML")
  sm.pvals <- anova(model.smoke.1, model.smoke.2)$p[2]
  print(sm.pvals)
  
  if(i==4){
    y.label="EP300 Expression (log CPM)"
    y.lim=c(7,8.5)
    y.text=8.5
  }
  graphTitle=paste(merrick.title,colnames(gsva.merrick.ep300)[i], sep=":")
  
  boxplot(gsva.vals~grade,main=graphTitle,
          col=c("#A2CD5A", "#E0EEC7", "#FFFFFF", "#F9BDBD","#F37C7C","#EE3B3B"),
          ylab=y.label, xlab="Grade", ylim=y.lim)
  text(x=4.5, y=y.text,label=paste("Anova P value =", formatC(grade.pval,format="e",digits=2), sep=" "))
  
  boxplot(gsva.vals~smoke, main=graphTitle, col=c("blue", "orange"),
          ylab=y.label, xlab="Smoking Status", ylim=y.lim)
  text(x=2.2, y=y.text,label=paste("Anova P value =", formatC(sm.pvals,format="e",digits=2), sep=" "))
}
dev.off()
