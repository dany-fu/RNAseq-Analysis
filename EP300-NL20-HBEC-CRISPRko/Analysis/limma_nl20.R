##################################################################################
# !! SOURCE LIMMA_TEST.R AND GSEA_FUNCTION.R BEFORE RUNNING CODE IN THIS FILE !! #
##################################################################################

# All NL20 samples
nl20.raw_counts = ep300.raw_counts[,13:24]
nl20.annot = rna.all.annot[rna.all.annot$CellLine=="NL20",]
dge.nl20 <- DGEList(counts=ep300.raw_counts[,13:24], samples=nl20.annot)
dim(dge.nl20) # [1] 60683    12

#################
# PREPROCESSING #
#################
dge.nl20.keep <- removeGenes(dge.nl20, dge.nl20$samples$Sample)
dim(dge.nl20.keep)  #14552
which(rownames(dge.nl20.keep) == "ENSG00000100393") 
dge.nl20.keep <- calcNormFactors(dge.nl20.keep, method = "TMM")
dge.nl20.keep.tmm <- cpm(dge.nl20.keep, log=TRUE) 

##############################
# BOXPLOT OF EP300 KNOCKDOWN #
##############################
ep300.rowNum <- which(rownames(dge.nl20.keep.tmm) == "ENSG00000100393") # 1753
nl20.bplot.data <- data.frame(Count.norm=unname(dge.nl20.keep.tmm[ep300.rowNum,]),
                              SampleID=names(dge.nl20.keep.tmm[ep300.rowNum,]),
                              Sample=nl20.annot$Sample)

tp53.ensId<-"ENSG00000141510"
notch1.ensId <- "ENSG00000148400"
tp53.row<-which(rownames(dge.nl20.keep.tmm) == tp53.ensId) # 5864
notch1.row<-which(rownames(dge.nl20.keep.tmm) == notch1.ensId) # 6513
nl20.tp53.bplot.data <- data.frame(Count.norm=unname(dge.nl20.keep.tmm[tp53.row,]),
                                   SampleID=names(dge.nl20.keep.tmm[tp53.row,]),
                                   Sample=nl20.annot$Sample)
nl20.notch1.bplot.data <- data.frame(Count.norm=unname(dge.nl20.keep.tmm[notch1.row,]),
                                   SampleID=names(dge.nl20.keep.tmm[notch1.row,]),
                                   Sample=nl20.annot$Sample)

pdf(file="Figures/boxplot-nl20.pdf")
par(mfrow=c(1,2))
print(ggboxplot(nl20.bplot.data, x="Sample", y="Count.norm", palette = "BuGn",fill="Sample",
          ylim = c(4.5, 7.6), ylab="Log2 EP300 Gene Level Counts per Million")+
        stat_compare_means(label = "p.signif",method = "t.test", ref.group = "NL20-CONTROL")+
        theme(legend.position = "none"))
print(ggboxplot(nl20.tp53.bplot.data, x="Sample", y="Count.norm", palette = "Reds", fill="Sample", 
                ylim = c(4.5, 7.6), ylab="Log2 TP53 Gene Level Counts per Million")+
        stat_compare_means(label = "p.signif",method = "t.test", ref.group = "NL20-CONTROL")+
        theme(legend.position = "none"))
dev.off()

############################
# MULTIDIMENSIONAL SCALING #
############################
makeMDSPlot(dge.nl20.keep.tmm, dge.nl20.keep)

#########
# LIMMA #
#########
nl20.fit_bayes <- transformAndFit(dge.nl20.keep) #voom will do cpm
summary(decideTests(nl20.fit_bayes))
# (Intercept) PerturbationEP300
# Down            29              1721
# NotSig        2933             11776
# Up           11590              1055
nl20.perturbation.table <- topTable(nl20.fit_bayes,coef=2, n=Inf,sort.by="p",adjust.method="BH", p.value=0.01,lfc=1.5) #560
which(rownames(nl20.perturbation.table) == "ENSG00000100393") # 546
dim(nl20.perturbation.table)
saveRDS(nl20.perturbation.table,file="nl20-degs.rds")

which(rownames(nl20.perturbation.table) == tp53.ensId) # 0

###########
# HEATMAP #
###########
geneNames <- convertEnsIdToGeneName(rownames(dge.nl20.keep.tmm))
i <- which(rownames(dge.nl20.keep.tmm) %in% rownames(nl20.perturbation.table))

nl20.ep300Genes.down <- rownames(subset(nl20.perturbation.table, logFC <0)) #479
nl20.ep300Genes.up <- rownames(subset(nl20.perturbation.table, logFC >0)) #81
#ep300Genes.down.names <- convertEnsIdToGeneName(ep300Genes.down)
#which(ep300Genes.down.names$hgnc_symbol == "EP300") # 469
#ep300Genes.up.names <- convertEnsIdToGeneName(ep300Genes.up)

colnames(dge.nl20.keep.tmm) <- c("NL20-CONTROL-1", "NL20-CONTROL-2", "NL20-CONTROL-3",
                                 "NL20-EP300-C11-1", "NL20-EP300-C11-2", "NL20-EP300-C11-3",
                                 "NL20-EP300-C18-1","NL20-EP300-C18-2","NL20-EP300-C18-3",
                                 "NL20-EP300-C20-1","NL20-EP300-C20-2","NL20-EP300-C20-3")

# Add sample annotations on the columns
annotation_col <- data.frame(as.factor(nl20.annot$Sample))
rownames(annotation_col) <- colnames(dge.nl20.keep.tmm)
colnames(annotation_col) <- 'Samples'

# Add gene annotations on the rows
ann_row <- data.frame(ep300.signature=c(rep("Upregulated", length(ep300Genes.up)), rep("Downregulated", length(ep300Genes.down))))
rownames(ann_row) <- c(ep300Genes.up, ep300Genes.down)
ann_row <- ann_row[match(rownames(dge.nl20.keep.tmm), rownames(ann_row)),,drop=FALSE]

# Add the label for the EP300 gene
geneLabels<- rep("", nrow(geneNames))
geneLabels[which(geneNames$hgnc_symbol == "EP300")] = "EP300"


breaksList = seq(-3, 3, by = .25)
pheatmap(dge.nl20.keep.tmm[i,], labels_row = geneLabels[i],
         # If you choose to not use the manual row scaling from above, select scale = "row"
         cluster_rows = T, cluster_cols = T, scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         annotation_row = ann_row,
         annotation_col = annotation_col,
         fontsize = 7,
         clustering_method = "ward.D2")

########
# GSEA #
########
nl20.pertTable.merge.geneNames <- merge(nl20.perturbation.table, geneNames, 
                                        by.x="row.names", by.y="ensembl_gene_id", sort=FALSE)
nl20.gsea.geneList = nl20.pertTable.merge.geneNames$t
names(nl20.gsea.geneList) = nl20.pertTable.merge.geneNames$hgnc_symbol
nl20.gsea.geneList <- nl20.gsea.geneList[!(is.na(names(nl20.gsea.geneList)) | names(nl20.gsea.geneList)=="")]
nl20.gsea.geneList = sort(nl20.gsea.geneList, decreasing = TRUE)
nl20.gsea.geneList = nl20.gsea.geneList[!duplicated(names(nl20.gsea.geneList))]

hallmark.geneSet <- "GSEA-genesets/h.all.v7.5.1.symbols.gmt"
oncogenic.geneSet <-"GSEA-genesets/c6.all.v7.5.1.symbols.gmt"
bio_processes.geneSet <- "GSEA-genesets/c5.go.bp.v7.5.1.symbols.gmt"
mol_func.geneSet <-"GSEA-genesets/c5.go.mf.v7.5.1.symbols.gmt"
nsclc.geneSet <- "GSEA-genesets/KEGG_NON_SMALL_CELL_LUNG_CANCER.v7.5.1.gmt"
kegg.geneSet <- "GSEA-genesets/c2.cp.kegg.v7.5.1.symbols.gmt"
GSEA(nl20.gsea.geneList, hallmark.geneSet, 0.05)




#################################
# ANALYSIS OF INDIVIDUAL CLONES #
#################################

#NL20 c11
nl20.c11.raw_counts = ep300.raw_counts[,13:18]
nl20.c11.annot = rna.all.annot[rna.all.annot$CellLine=="NL20" & (rna.all.annot$Sample=="NL20-EP300-11" | rna.all.annot$Sample=="NL20-CONTROL"),]
dge.nl20.c11.withAnnot <- DGEList(counts=nl20.c11.raw_counts, samples=nl20.c11.annot)
dim(dge.nl20.c11.withAnnot) # [1] 60683    6
dge.nl20.c11.withAnnot.keep <- removeGenes(dge.nl20.c11.withAnnot,dge.nl20.c11.withAnnot$samples$Sample)
dim(dge.nl20.c11.withAnnot.keep) # 13818
which(rownames(dge.nl20.c11.withAnnot.keep) == "ENSG00000100393") 
dge.nl20.c11.withAnnot.keep <- calcNormFactors(dge.nl20.c11.withAnnot.keep, method = "TMM")
nl20.c11.fit_bayes <- transformAndFit(dge.nl20.c11.withAnnot.keep)
summary(decideTests(nl20.c11.fit_bayes))
# (Intercept) PerturbationEP300
# Down           361              4745
# NotSig        1139              4338
# Up           12318              4735
nl20.c11.perturbation.table <- topTable(nl20.c11.fit_bayes,coef=2, n=Inf,sort.by="p",adjust.method="BH", p.value=0.01, lfc=1.5) #1860
which(rownames(nl20.c11.perturbation.table) == "ENSG00000100393") # 281
which(rownames(nl20.c11.perturbation.table) == tp53.ensId) #0
dim(nl20.c11.perturbation.table)

nl20.c11.pertTable.merge.geneNames <- merge(nl20.c11.perturbation.table, geneNames, 
                                            by.x="row.names", by.y="ensembl_gene_id", sort=FALSE)
nl20.c11.gsea.geneList = nl20.c11.pertTable.merge.geneNames$t
names(nl20.c11.gsea.geneList) = nl20.c11.pertTable.merge.geneNames$hgnc_symbol
nl20.c11.gsea.geneList <- nl20.c11.gsea.geneList[!(is.na(names(nl20.c11.gsea.geneList)) | names(nl20.c11.gsea.geneList)=="")]
nl20.c11.gsea.geneList = sort(nl20.c11.gsea.geneList, decreasing = TRUE)
nl20.c11.gsea.geneList = nl20.c11.gsea.geneList[!duplicated(names(nl20.c11.gsea.geneList))]
GSEA(nl20.c11.gsea.geneList, hallmark.geneSet, 0.05)

#NL20 c18
nl20.c18.raw_counts = ep300.raw_counts[,c(13:15, 19:21)]
nl20.c18.annot = rna.all.annot[rna.all.annot$CellLine=="NL20" & (rna.all.annot$Sample=="NL20-EP300-18" | rna.all.annot$Sample=="NL20-CONTROL"),]
dge.nl20.c18.withAnnot <- DGEList(counts=nl20.c18.raw_counts, samples=nl20.c18.annot)
dim(dge.nl20.c18.withAnnot) # [1] 60683    6
dge.nl20.c18.withAnnot.keep <- removeGenes(dge.nl20.c18.withAnnot, dge.nl20.c18.withAnnot$samples$Sample)
dim(dge.nl20.c18.withAnnot.keep) #13858
which(rownames(dge.nl20.c18.withAnnot.keep) == "ENSG00000100393") 
dge.nl20.c18.withAnnot.keep <- calcNormFactors(dge.nl20.c18.withAnnot.keep, method = "TMM")
nl20.c18.fit_bayes <- transformAndFit(dge.nl20.c18.withAnnot.keep)
summary(decideTests(nl20.c18.fit_bayes))
# (Intercept) PerturbationEP300
# Down           391              4858
# NotSig        1099              4238
# Up           12368              4762
nl20.c18.perturbation.table <- topTable(nl20.c18.fit_bayes,coef=2, n=Inf,sort.by="p",adjust.method="BH", p.value=0.01, lfc=1.5) #1869
which(rownames(nl20.c18.perturbation.table) == "ENSG00000100393") #174
which(rownames(nl20.c18.perturbation.table) == tp53.ensId) #0
dim(nl20.c18.perturbation.table)

nl20.c18.ep300Genes.down <- rownames(subset(nl20.c18.perturbation.table, logFC <0)) #1065
nl20.c18.ep300Genes.up <- rownames(subset(nl20.c18.perturbation.table, logFC >0))  #804

nl20.c18.pertTable.merge.geneNames <- merge(nl20.c18.perturbation.table, geneNames, 
                                            by.x="row.names", by.y="ensembl_gene_id", sort=FALSE)
nl20.c18.gsea.geneList = nl20.c18.pertTable.merge.geneNames$t
names(nl20.c18.gsea.geneList) = nl20.c18.pertTable.merge.geneNames$hgnc_symbol
nl20.c18.gsea.geneList <- nl20.c18.gsea.geneList[!(is.na(names(nl20.c18.gsea.geneList)) | names(nl20.c18.gsea.geneList)=="")]
nl20.c18.gsea.geneList = sort(nl20.c18.gsea.geneList, decreasing = TRUE)
nl20.c18.gsea.geneList = nl20.c18.gsea.geneList[!duplicated(names(nl20.c18.gsea.geneList))]
GSEA(nl20.c18.gsea.geneList, hallmark.geneSet, 0.05)

hallmark.geneSet <- "GSEA-genesets/h.all.v7.5.1.symbols.gmt"
oncogenic.geneSet <-"GSEA-genesets/c6.all.v7.5.1.symbols.gmt"
bio_processes.geneSet <- "GSEA-genesets/c5.go.bp.v7.5.1.symbols.gmt"
mol_func.geneSet <-"GSEA-genesets/c5.go.mf.v7.5.1.symbols.gmt"
nsclc.geneSet <- "GSEA-genesets/KEGG_NON_SMALL_CELL_LUNG_CANCER.v7.5.1.gmt"
kegg.geneSet <- "GSEA-genesets/c2.cp.kegg.v7.5.1.symbols.gmt"

#NL20 c20
nl20.c20.raw_counts = ep300.raw_counts[,c(13:15, 22:24)]
nl20.c20.annot = rna.all.annot[rna.all.annot$CellLine=="NL20" & (rna.all.annot$Sample=="NL20-EP300-20" | rna.all.annot$Sample=="NL20-CONTROL"),]
dge.nl20.c20.withAnnot <- DGEList(counts=nl20.c20.raw_counts, samples=nl20.c20.annot)
dim(dge.nl20.c20.withAnnot) # [1] 60683    6
dge.nl20.c20.withAnnot.keep <- removeGenes(dge.nl20.c20.withAnnot, dge.nl20.c20.withAnnot$samples$Sample)
dim(dge.nl20.c20.withAnnot.keep) #13343
which(rownames(dge.nl20.c20.withAnnot.keep) == "ENSG00000100393") 
dge.nl20.c20.withAnnot.keep <- calcNormFactors(dge.nl20.c20.withAnnot.keep, method = "TMM")
nl20.c20.fit_bayes <- transformAndFit(dge.nl20.c20.withAnnot.keep)
summary(decideTests(nl20.c20.fit_bayes))
# (Intercept) PerturbationEP300
# Down           156              3846
# NotSig         917              5595
# Up           12270              3902
nl20.c20.perturbation.table <- topTable(nl20.c20.fit_bayes,coef=2, n=Inf,sort.by="p",adjust.method="BH", p.value=0.01, lfc=1.5) #0.5 is the greatest fc
which(rownames(nl20.c20.perturbation.table) == "ENSG00000100393") #1689
dim(nl20.c20.perturbation.table)
which(rownames(nl20.c20.perturbation.table) == tp53.ensId) #0







###################
# Just c11 and 18 #
###################
nl20.best.raw_counts = ep300.raw_counts[,c(13:21)]
nl20.best.annot = rna.all.annot[rna.all.annot$CellLine=="NL20" & (rna.all.annot$Sample!="NL20-EP300-20"),]
dge.nl20.best <- DGEList(counts=nl20.best.raw_counts, samples=nl20.best.annot)
dim(dge.nl20.best) # [1] 60683    9
dge.nl20.best.keep <- removeGenes(dge.nl20.best, dge.nl20.best$samples$Sample)
dim(dge.nl20.best.keep) #14350
which(rownames(dge.nl20.best.keep) == "ENSG00000100393") #1739
dge.nl20.best.keep <- calcNormFactors(dge.nl20.best.keep, method = "TMM")
nl20.best.fit_bayes <- transformAndFit(dge.nl20.best.keep)
summary(decideTests(nl20.best.fit_bayes))
# (Intercept) PerturbationEP300
# Down           169              3310
# NotSig        2265              8419
# Up           11916              2621
nl20.best.perturbation.table <- topTable(nl20.best.fit_bayes,coef=2, n=Inf,sort.by="p",adjust.method="BH", p.value=0.01, lfc=1.5) #1191
which(rownames(nl20.best.perturbation.table) == "ENSG00000100393") #444

saveRDS(nl20.best.perturbation.table,file="nl20-best-degs.rds")


###########
# HEATMAP #
###########
nl20.best.lcpm <- cpm(dge.nl20.best, log=TRUE)
dim(nl20.best.lcpm)
geneNames <- convertEnsIdToGeneName(rownames(nl20.best.lcpm))
i <- which(rownames(nl20.best.lcpm) %in% rownames(nl20.best.perturbation.table))

nl20.best.ep300Genes.down <- rownames(subset(nl20.best.perturbation.table, logFC <0)) #847
nl20.best.ep300Genes.up <- rownames(subset(nl20.best.perturbation.table, logFC >0)) #344
#ep300Genes.down.names <- convertEnsIdToGeneName(ep300Genes.down)
#which(ep300Genes.down.names$hgnc_symbol == "EP300") # 469
#ep300Genes.up.names <- convertEnsIdToGeneName(ep300Genes.up)

# Add sample annotations on the columns
annotation_col <- data.frame(as.factor(nl20.best.annot$Sample))
rownames(annotation_col) <- colnames(nl20.best.lcpm)
colnames(annotation_col) <- 'Samples'

# Add gene annotations on the rows
ann_row <- data.frame(ep300.signature=c(rep("Upregulated", length(ep300Genes.up)), rep("Downregulated", length(ep300Genes.down))))
rownames(ann_row) <- c(ep300Genes.up, ep300Genes.down)
ann_row <- ann_row[match(rownames(nl20.best.lcpm), rownames(ann_row)),,drop=FALSE]

# Add the label for the EP300 gene
which(geneNames$hgnc_symbol == "EP300") # 2294
geneLabels<- rep("", 60683)
geneLabels[2294] = "EP300"

breaksList = seq(-3, 3, by = .25)
pheatmap(nl20.best.lcpm[i,], labels_row = geneLabels[i],
         # If you choose to not use the manual row scaling from above, select scale = "row"
         cluster_rows = T, cluster_cols = T, scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         annotation_row = ann_row,
         annotation_col = annotation_col,
         fontsize = 7,
         clustering_method = "ward.D2")

########
# GSEA #
########
nl20.best.perturbation.table<-readRDS("GSVA/nl20-best-degs.rds")
nl20.best.pertTable.merge.geneNames <- merge(nl20.best.perturbation.table, geneNames, 
                                             by.x="row.names", by.y="ensembl_gene_id", sort=FALSE)
nl20.best.gsea.geneList = nl20.best.pertTable.merge.geneNames$t
names(nl20.best.gsea.geneList) = nl20.best.pertTable.merge.geneNames$hgnc_symbol
nl20.best.gsea.geneList <- nl20.best.gsea.geneList[!(is.na(names(nl20.best.gsea.geneList)) | names(nl20.best.gsea.geneList)=="")]
nl20.best.gsea.geneList = sort(nl20.best.gsea.geneList, decreasing = TRUE)
nl20.best.gsea.geneList = nl20.best.gsea.geneList[!duplicated(names(nl20.best.gsea.geneList))]
GSEA(nl20.best.gsea.geneList, hallmark.geneSet, 0.05)

###############################
# Make geneset from c11 and 18#
###############################
nl20.best.gsea.geneList.up <- subset(nl20.best.gsea.geneList, nl20.best.gsea.geneList>0)
df.geneset.up<- data.frame(names(nl20.best.gsea.geneList.up))
write.table(x = t(df.geneset.up),
            file = "nl20-geneset.gmt",
            row.names = F,
            col.names = F,
            quote = F,
            sep = '\t')

nl20.best.gsea.geneList.down <- subset(nl20.best.gsea.geneList, nl20.best.gsea.geneList<0)
df.geneset.down<- data.frame(names(nl20.best.gsea.geneList.down))
write.table(x = t(df.geneset.down),
            file = "nl20-geneset.gmt",
            row.names = F,
            col.names = F,
            quote = F,
            sep = '\t',
            append=T)
