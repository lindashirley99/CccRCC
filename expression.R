###This R script is used to perform gene expression analysis
###code starts########

#load in data
#######################################################################################
expr= read.csv("RNAseq_expr_RSEM.csv", h=T, stringsAsFactors = F, row.names = 1, check.names = F)
expr = expr[,order(colnames(expr))]
samples = c(grep("N$",colnames(expr),value = T,invert = T),grep("N$",colnames(expr),value = T))
expr = expr[,samples]
anno = read.table("Sample_info-2017-12-14.txt", h=T, sep = "\t",stringsAsFactors = F, row.names = 1, check.names = F)
cc = row.names(anno[anno[,2]=="Clear_cell",])
anno = anno[cc,]
anno[,"Race"] = rep("asian",nrow(anno))
anno[anno[,"Stage"]%in% c("T1a","T1b","T2a","T2b"),"Stage_group"] = "T1T2"
anno[anno[,"Stage"]%in% c("T3","T3a","T3b","T3c","T4"),"Stage_group"] = "T3T4"
anno[anno[,"Grade"]%in% c("Grade1","Grade2"),"Grade_group"] = "G1G2"
anno[anno[,"Grade"]%in% c("Grade3","Grade4"),"Grade_group"] = "G3G4"

expr = expr[,cc]
write.table(expr, file = "RNAseq_expr_RSEM_CccRCC.txt", sep = "\t")

#TCGA ccRCC gene expression data
TCGA_T = read.table("./TCGA/TCGA_T_gene.KIRC.info", h=T, stringsAsFactors = F, 
                    sep="\t", row.names = 1, check.names = F)
TCGA_T[,1] = NULL
TCGA_N = read.table("./TCGA/TCGA_N_gene.KIRC.info", h=T, stringsAsFactors = F, 
                    sep="\t", row.names = 1, check.names = F)
TCGA_N[,1] = NULL
genes = row.names(expr)
TCGA_T = TCGA_T[genes,]
TCGA_N = TCGA_N[genes,]
anno_TCGA = read.table("./TCGA/TCGA_kidney_clinical-2017-11-29.txt", h=T, 
                     stringsAsFactors = F, sep="\t", check.names = F)
anno_TCGA = anno_TCGA[anno_TCGA[,3] == "Kidney Clear Cell Renal Carcinoma",]
anno_TCGA[,2] = tolower(anno_TCGA[,2])
clinical =  read.table("./TCGA/TCGA_KIRC_clinical_2018-07-05.tsv", h=T, 
                       stringsAsFactors = F, sep="\t", check.names = F)
anno_TCGA = merge(anno_TCGA[,c(1:6,10:12)], clinical[,c(1,2,14,16,25)],by.x=2,by.y=1)
anno_TCGA[anno_TCGA[,"stage"]%in% c("StageI","StageII"),"Stage_group"] = "T1T2"
anno_TCGA[anno_TCGA[,"stage"]%in% c("StageIII","StageIV"),"Stage_group"] = "T3T4"
anno_TCGA[anno_TCGA[,"grade"]%in% c("G1","G2"),"Grade_group"] = "G1G2"
anno_TCGA[anno_TCGA[,"grade"]%in% c("G3","G4"),"Grade_group"] = "G3G4"
anno_TCGA[anno_TCGA[,"grade"]%in% c("GX"),"Grade_group"] = NA
anno_TCGA[anno_TCGA[,"metastasis"]%in% c("M0"),"Metastatic"] = "N"
anno_TCGA[anno_TCGA[,"metastasis"]%in% c("M1","MX"),"Metastatic"] = "Y"
anno_TCGA[anno_TCGA[,"metastasis"]=="","Metastatic"] = NA
anno_TCGA[anno_TCGA[,"histology"]=="Kidney Clear Cell Renal Carcinoma","Type"] = "Clear_cell"
anno_TCGA[anno_TCGA[,"histology"]=="Kidney Papillary Renal Cell Carcinoma","Type"] = "Pappilary"
anno_TCGA[anno_TCGA[,"histology"]=="Kidney Chromophobe","Type"] = "Chromophobe"
anno_TCGA[anno_TCGA[,"gender"]=="FEMALE","gender"] = "F"
anno_TCGA[anno_TCGA[,"gender"]=="MALE","gender"] = "M"
row.names(anno_TCGA)=anno_TCGA[,2]
anno_TCGA[,c(2,3,10)]= NULL
anno_TCGA[,c(2,4:8,11:14)] <- lapply(anno_TCGA[,c(2,4:8,11:14)], as.factor)

samples_T=colnames(TCGA_T)[colnames(TCGA_T) %in% row.names(anno_TCGA)]
samples_N=colnames(TCGA_N)[colnames(TCGA_N) %in% row.names(anno_TCGA)]
TCGA_T=TCGA_T[,samples_T] 
TCGA_N=TCGA_N[,samples_N]

#Table 1: clinical statistics
library("plyr")
count(anno_TCGA[samples_T,8]) 
mean(anno_TCGA[samples_T,3]) 
############################################################################


temp1=cbind(rep("Tumor", 533),anno_TCGA[colnames(TCGA_T),c(14,2:5,11,6,12,7)])
colnames(temp1) = colnames(anno)[c(1:4,22,5,23,6,24,8)]
temp2=cbind(rep("Normal", 72),anno_TCGA[colnames(TCGA_N),c(14,2:5,11,6,12,7)])
colnames(temp2) = colnames(anno)[c(1:4,22,5,23,6,24,8)]
pheno = rbind(anno[colnames(expr),c(1:4,22,5,23,6,24,8)],temp1, temp2)
pheno[,"Dataset"] = c(rep("C",66),rep("T",605))
pheno[,"Race"] = tolower(pheno[,"Race"])
pheno[pheno[,"Race"]=="","Race"] = "NA"
pheno[pheno[,"Race"]=="black or african american","Race"] = "B"
pheno[pheno[,"Race"]=="white","Race"] = "W"
pheno[pheno[,"Race"]=="asian","Race"] = "A"


dat = cbind(expr,TCGA_T,TCGA_N)
samples = row.names(pheno)
colnames(dat) = samples

#remove batch effect
##############################################################################
library("COMBAT")
library("sva")
mod = model.matrix(~as.factor(Group), data=pheno[,1:5])
mod0 = model.matrix(~1,data=pheno[,1:5])
n.sv = num.sv(dat,mod,method="leek") #estimate batch and other artifacts
svobj = sva(as.matrix(dat),mod,mod0,n.sv=n.sv)

combat_dat = ComBat(as.matrix(dat), pheno$Dataset, mod=mod, par.prior = TRUE, prior.plots = FALSE)
boxplot(log2(dat[,50:99]))
boxplot(log2(combat_dat[,50:99]))
#############################################################################
#Fig.1b PCA of ccRCC for CKC and TCGA colored by race

dat = cbind(expr[,colnames(expr_t)],TCGA_T)
samples = colnames(dat)
sample_anno = pheno[samples,]

library("Rtsne")
tsne <- Rtsne(as.matrix(t(dat)), dims = 2, perplexity=30, verbose=TRUE, max_iter = 500,check_duplicates=F)
pdf("tsne CKC TCGA with text labels.pdf",20,20)
colors = rainbow(length(unique(sample_anno$Race)))
names(colors) = unique(sample_anno$Race)
plot(tsne$Y, t='n', main="tsne of ccRCC patients in CKC and TCGA")
text(tsne$Y, labels=sample_anno$Dataset, col=colors[sample_anno$Race],cex = 1.5)  
legend("topright", legend=names(colors), col=unique(colors), cex=1.8, pch = 19)
dev.off()

pdf("tsne of ccRCC patients in CKC and TCGA KIRC.pdf",20,20)
colors = rainbow(length(unique(sample_anno$Race)))
names(colors) = unique(sample_anno$Race)
plot(tsne$Y, t='n', main="tsne of ccRCC patients in CKC and TCGA KIRC")
points(tsne$Y, pch = c(rep(21,55),rep(24,533)), col=colors[sample_anno$Race],cex = 1.5) 
legend("topright", legend=names(colors), col=unique(colors), cex=1.8, pch = 19)
legend("bottomright", legend=c("CKC","TCGA"), pch = c(21,24), cex=1.8)
dev.off()

library(ggplot2)
dat=expr_t[rowSums(expr_t) > 0,]
pca <- prcomp(as.matrix(t(dat)), center = TRUE, scale. = TRUE)
samples=colnames(expr_t)
eigs <- pca$sdev^2
proVar = eigs/sum(eigs)*100
scores = as.data.frame(pca$x)
scores$Stage_group=anno[samples,"Stage_group"]
scores$Grade_group=anno[samples,"Grade_group"]
pdf("PCA_CccRCC_grade.pdf",10,8)
p <- ggplot(data = scores, aes(x = PC1, y = PC2))+ labs(x = paste("PC1(",format(proVar[1], digits=2, nsmall=2),"%)",sep = ""), y = paste("PC2(",format(proVar[2], digits=2, nsmall=2),"%)",sep = ""))
p <- p + geom_point(aes(x = PC1, y = PC2,fill=Grade_group),color="black",shape=21,size=5, stat = "identity") 
p + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("PCA_CccRCC_stage.pdf",10,8)
p <- ggplot(data = scores, aes(x = PC1, y = PC2))+ labs(x = paste("PC1(",format(proVar[1], digits=2, nsmall=2),"%)",sep = ""), y = paste("PC2(",format(proVar[2], digits=2, nsmall=2),"%)",sep = ""))
p <- p + geom_point(aes(x = PC1, y = PC2,fill=Stage_group),color="black",shape=21,size=5, stat = "identity") 
p + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
##########################################################################################

expr_g = expr
expr_g[,"Gene"] = sub("\\|.+", "", row.names(expr_g))
expr_g = expr_g[!duplicated(expr_g[,"Gene"]),]
row.names(expr_g) = expr_g[,"Gene"]
expr_g[,"Gene"] = NULL

TCGA_T_g = TCGA_T
TCGA_T_g[,"Gene"] = sub("\\|.+", "", row.names(TCGA_T_g))
TCGA_T_g = TCGA_T_g[!duplicated(TCGA_T_g[,"Gene"]),]
row.names(TCGA_T_g) = TCGA_T_g[,"Gene"]
TCGA_T_g[,"Gene"] = NULL


#DE genes asian vs. white in ccRCC
######################################################################
sample_anno_cc = subset(pheno, Group == "Tumor" & Type == "Clear_cell")
sample_anno = subset(sample_anno_cc,Race == "A" | Race == "W")
samples = row.names(sample_anno)
dat = cbind(expr_g,TCGA_T_g)
dat = as.matrix(dat[,samples])
library("matrixStats")
dat_ex = dat[rowQuantiles(as.matrix(dat),probs = 0.2) > 100,]
library("limma")
group <- factor(as.vector(sample_anno[,5]))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(log2(dat_ex+1), design)
contrast.matrix <- makeContrasts("A-W",levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2,number=nrow(fit2), adjust="BH",coef="A-W")
PvalueCutoff <- 0.05
logFCcutoff <- log2(2)
de <- topTable(fit2, number=Inf, coef="A-W",p.value=PvalueCutoff,adjust="BH",lfc = logFCcutoff )
pdf("Heatmap-DE genes ccRCC asian vs. white.pdf",14,12)
library("pheatmap")
par(mar=c(5.1, 5.1, 5.1, 9.1), xpd=TRUE)
pheatmap(dat[row.names(de),],cluster_rows=T,cluster_cols=T,treeheight_col = 15,scale = "row", annotation_col = sample_anno,
         annotation_legend = T,legend=T, show_rownames = T,show_colnames = F,
         fontsize=5,  color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

#GSVA
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(RColorBrewer)
require("biomaRt")
data(c2BroadSets)
canonicalC2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                                      + grep("^REACTOME", names(c2BroadSets)),
                                      + grep("^BIOCARTA", names(c2BroadSets)))]
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
gs=rownames(dat)
gsInfo=getBM(attributes=c('hgnc_symbol','entrezgene','gene_biotype'),filters = 'hgnc_symbol', 
             values = gs,mart=ensembl)
eids=gsInfo[[2]]
names(eids)=gsInfo[[1]]
eids=eids[gs]
eids=eids[!is.na(eids)]
rma_nna=as.matrix(dat[names(eids),])
rownames(rma_nna)=eids
temp = gsva(rma_nna,canonicalC2BroadSets,min.sz=5,max.sz=500,rnaseq=T,mx.diff=TRUE,method="gsva",
            annotation="org.Hs.eg.db")
expr_gs = temp$es.obs
fit <- lmFit(expr_gs, design)
contrast.matrix <- makeContrasts("A-W", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
allGeneSets <- topTable(fit2,number=nrow(fit2), adjust="BH",coef="A-W")
PvalueCutoff <- 0.05
logFCcutoff <- log2(1.3)
DEgeneSets <- topTable(fit2, coef="A-W", number=Inf, p.value=PvalueCutoff,adjust="BH",lfc = logFCcutoff)
DEgeneSets
GSVAsco <- expr_gs[rownames(DEgeneSets), ]
pdf("GSVA-DEGeneSets-ccRCC asian vs. white.pdf", 16,10)
pheatmap(GSVAsco,cluster_rows=T,cluster_cols=T, treeheight_row=15,treeheight_col = 15, annotation_col= sample_anno,
         scale = "row", annotation_legend=T,legend=T,show_rownames = T,show_colnames = F,
         cellheight=10, cellwidth=1,fontsize=8, fontsize_row = 8, 
         main = "Adjusted P-value < 0.05;|FC|>1.3",
         color = colorRampPalette(c("white", "yellow", "red"))(100))
dev.off()

temp = gsva(rma_nna,canonicalC2BroadSets,min.sz=5,max.sz=500,rnaseq=T,mx.diff=TRUE,method="ssgsea",
            annotation="org.Hs.eg.db")

pdf("ssGSEA-ccRCC asian vs. white.pdf", 16,10)
library("pheatmap")
pheatmap(temp,cluster_rows=T,cluster_cols=T, treeheight_row=15,treeheight_col = 15, annotation_col= sample_anno,
         scale = "row", annotation_legend=T,legend=T,show_rownames = T,show_colnames = F,
         cellheight=10, cellwidth=1,fontsize=8, fontsize_row = 8, 
         main = "ssGSEA",color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()
eids=gsInfo[[2]]
names(eids)=gsInfo[[1]]
eids=eids[gs]
eids=eids[!is.na(eids)]
rma_nna=as.matrix(expr_g[names(eids),])
rownames(rma_nna)=eids
#######################################################################################

#detect DE genes/pathways 55 tumor vs. 11 normal
#############################################################
library("limma")
dat = expr_g
group <- factor(as.vector(anno[,"Group"]))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(log2(dat+1), design)
contrast.matrix <- makeContrasts("Tumor-Normal", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2,number=nrow(fit2), adjust="BH",coef="Tumor-Normal")
PvalueCutoff <- 0.0001
logFCcutoff <- log2(4)
de <- topTable(fit2, coef="Tumor-Normal", number=Inf, p.value=PvalueCutoff,adjust="BH",lfc = logFCcutoff )
pdf("Heatmap-DE genes T vs N.pdf",14,12)
library("pheatmap")
par(mar=c(5.1, 5.1, 5.1, 9.1), xpd=TRUE)
pheatmap(dat[row.names(de),],cluster_rows=T,cluster_cols=T,treeheight_col = 15,scale = "row", annotation_col = anno[,1:2],
         annotation_legend = T,legend=T, show_rownames = F,show_colnames = T,
         fontsize=8,  color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

#GSVA
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(RColorBrewer)
require("biomaRt")
data(c2BroadSets)
canonicalC2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                                      + grep("^REACTOME", names(c2BroadSets)),
                                      + grep("^BIOCARTA", names(c2BroadSets)))]
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
gs=rownames(dat)
gsInfo=getBM(attributes=c('hgnc_symbol','entrezgene_id','gene_biotype'),filters = 'hgnc_symbol', 
             values = gs,mart=ensembl)
eids=gsInfo[[2]]
names(eids)=gsInfo[[1]]
eids=eids[gs]
eids=eids[!is.na(eids)]
rma_nna=as.matrix(dat[names(eids),])
rownames(rma_nna)=eids
temp = gsva(rma_nna,canonicalC2BroadSets,min.sz=5,max.sz=500,rnaseq=T,mx.diff=TRUE,method="gsva",annotation="org.Hs.eg.db")
expr_gs = temp$es.obs
group <- factor(as.vector(anno[,"Group"]))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_gs, design)
contrast.matrix <- makeContrasts("Tumor-Normal", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
allGeneSets <- topTable(fit2,number=nrow(fit2), adjust="BH",coef="Tumor-Normal")
PvalueCutoff <- 0.05
logFCcutoff <- log2(1.8)
DEgeneSets <- topTable(fit2, coef="Tumor-Normal", number=Inf, p.value=PvalueCutoff,adjust="BH",lfc = logFCcutoff )
DEgeneSets
GSVAsco <- expr_gs[rownames(DEgeneSets), ]
pdf("GSVA-DEGeneSets-T vs N.pdf", 18,12)
pheatmap(GSVAsco,cluster_rows=T,cluster_cols=T, treeheight_row=15,treeheight_col = 15, annotation_col= anno[,1,drop=F],
         scale = "row", annotation_legend=T,legend=T,show_rownames = T,show_colnames = T,
         cellheight=10, cellwidth=10,fontsize=8, fontsize_row = 8, 
         main = "Adjusted P-value < 0.05; |FC| > 1.6",
         color = colorRampPalette(c("white", "yellow", "red"))(100))
dev.off()

#DE genes/pathways across stages
library("limma")
expr_t= as.matrix(expr_g[,1:55])
samples = colnames(expr_t)
stages = anno[samples,c(5,6,23,24)]
stages = na.omit(stages)
samples = row.names(stages)
dat = expr_t[,samples]
library("matrixStats")
dat_ex = dat[rowQuantiles(as.matrix(dat),probs = 0.2) > 100,]
group <- factor(as.vector(stages[,"Stage_group"]))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(log2(dat_ex + 1), design)
contrast.matrix <- makeContrasts("T1T2-T3T4", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2,number=nrow(fit2), adjust="BH",coef="T1T2-T3T4")
PvalueCutoff <- 0.05
logFCcutoff <- log2(2)
de <- topTable(fit2, number=Inf, coef="T1T2-T3T4",p.value=PvalueCutoff,adjust="none",lfc = logFCcutoff )

pdf("Heatmap-DE genes T12 vs. T34.pdf",14,12)
library("pheatmap")
par(mar=c(5.1, 5.1, 5.1, 9.1), xpd=TRUE)
pheatmap(dat[row.names(de),],cluster_rows=T,cluster_cols=T,treeheight_col = 15,scale = "row", annotation_col = stages,
         annotation_legend = T,legend=T, show_rownames = T,show_colnames = T,
         fontsize=8,  color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

#GSVA
gs=rownames(dat)
gsInfo=getBM(attributes=c('hgnc_symbol','entrezgene_id','gene_biotype'),filters = 'hgnc_symbol', 
             values = gs,mart=ensembl)
eids=gsInfo[[2]]
names(eids)=gsInfo[[1]]
eids=eids[gs]
eids=eids[!is.na(eids)]
rma_nna=as.matrix(dat[names(eids),])
rownames(rma_nna)=eids
expr_gs = gsva(rma_nna,canonicalC2BroadSets,min.sz=5,max.sz=500,mx.diff=TRUE,method="gsva",annotation="org.Hs.eg.db")

group <- factor(as.vector(stages[,"Stage_group"]))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_gs, design)
contrast.matrix <- makeContrasts("T1T2-T3T4", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
allGeneSets <- topTable(fit2,number=nrow(fit2), adjust="BH",coef="T1T2-T3T4")
PvalueCutoff <- 0.05
logFCcutoff <- log2(1.3)
DEgeneSets <- topTable(fit2, coef="T1T2-T3T4", number=Inf, p.value=PvalueCutoff,adjust="none")
DEgeneSets
GSVAsco <- expr_gs[rownames(DEgeneSets), ]
pdf("GSVA-DEGeneSets-T12 vs. T34.pdf", 18,12)
library("pheatmap")
pheatmap(GSVAsco,cluster_rows=T,cluster_cols=T, treeheight_row=15,treeheight_col = 15, annotation_col= stages[,3,drop=F],
         scale = "row", annotation_legend=T,legend=T,show_rownames = T,show_colnames = T,
         cellheight=10, cellwidth=10,fontsize=8, fontsize_row = 8, 
         main = "|FC| > 1.3, P-value < 0.05",
         color = colorRampPalette(c("white", "yellow", "red"))(100))
dev.off()

#DE across grades
samples = colnames(expr_t)
grades = anno[samples,c(5,6,23,24)]
grades = na.omit(grades)
samples = row.names(grades)
dat = expr_t[,samples]
library("matrixStats")
dat_ex = dat[rowQuantiles(as.matrix(dat),probs = 0.2) > 100,]
group <- factor(as.vector(grades[,"Grade_group"]))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(log2(dat_ex + 1), design)
contrast.matrix <- makeContrasts("G1G2-G3G4", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2,number=nrow(fit2), adjust="BH",coef="G1G2-G3G4")
PvalueCutoff <- 0.05
logFCcutoff <- log2(2)
de <- topTable(fit2, number=Inf, coef="G1G2-G3G4",p.value=PvalueCutoff,adjust="BH",lfc = logFCcutoff )

pdf("Heatmap-DE genes G12 vs. G34.pdf",14,12)
library("pheatmap")
par(mar=c(5.1, 5.1, 5.1, 9.1), xpd=TRUE)
pheatmap(dat[row.names(de),],cluster_rows=T,cluster_cols=T,treeheight_col = 15,scale = "row", annotation_col = grades,
         annotation_legend = T,legend=T, show_rownames = T,show_colnames = T,
         fontsize=5,  color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

#GSVA
gs=rownames(dat)
gsInfo=getBM(attributes=c('hgnc_symbol','entrezgene_id','gene_biotype'),filters = 'hgnc_symbol', 
             values = gs,mart=ensembl)
eids=gsInfo[[2]]
names(eids)=gsInfo[[1]]
eids=eids[gs]
eids=eids[!is.na(eids)]
rma_nna=as.matrix(dat[names(eids),])
rownames(rma_nna) = eids
expr_gs = gsva(rma_nna,canonicalC2BroadSets,min.sz=5,max.sz=500,mx.diff=TRUE,method="gsva",annotation="org.Hs.eg.db")

group <- factor(as.vector(grades[,"Grade_group"]))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_gs, design)
contrast.matrix <- makeContrasts("G1G2-G3G4", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
allGeneSets <- topTable(fit2,number=nrow(fit2), adjust="BH",coef="G1G2-G3G4")
PvalueCutoff <- 0.05
logFCcutoff <- log2(1.3)
DEgeneSets <- topTable(fit2, coef="G1G2-G3G4", number=Inf, p.value=PvalueCutoff,adjust="BH",lfc = logFCcutoff)
DEgeneSets
GSVAsco <- expr_gs[rownames(DEgeneSets), ]
pdf("GSVA-DEGeneSets-G12 vs. G34.pdf", 18,12)
pheatmap(GSVAsco,cluster_rows=T,cluster_cols=T, treeheight_row=15,treeheight_col = 15, annotation_col= grades,
         scale = "row", annotation_legend=T,legend=T,show_rownames = T,show_colnames = T,
         cellheight=10, cellwidth=10,fontsize=8, fontsize_row = 8, 
         main = "Adjusted P-value < 0.05;|FC|>1.3",
         color = colorRampPalette(c("white", "yellow", "red"))(100))
dev.off()

#DE across hypertension or Diabetes
#############################################################
samples = colnames(expr_t)
meta = anno[samples,14:15]
meta = na.omit(meta)
samples = row.names(meta)
dat = expr_t[,samples]
library("matrixStats")
dat_ex = dat[rowQuantiles(as.matrix(dat),probs = 0.2) > 100,]
group <- factor(as.vector(meta[,"Drink"]))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(log2(dat_ex + 1), design)
contrast.matrix <- makeContrasts("Y-N", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2,number=nrow(fit2), adjust="BH",coef="Y-N")
PvalueCutoff <- 0.05
logFCcutoff <- log2(2)
de <- topTable(fit2, number=Inf, coef="Y-N",p.value=PvalueCutoff,adjust="none",lfc = logFCcutoff )

pdf("Heatmap-DE genes drinking vs. non-drinking tumors.pdf",14,12)
library("pheatmap")
par(mar=c(5.1, 5.1, 5.1, 9.1), xpd=TRUE)
pheatmap(dat[row.names(de),],cluster_rows=T,cluster_cols=T,treeheight_col = 15,scale = "row", annotation_col = meta,
         annotation_legend = T,legend=T, show_rownames = T,show_colnames = T,
         fontsize=5,  color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()


gs=rownames(dat)
gsInfo=getBM(attributes=c('hgnc_symbol','entrezgene','gene_biotype'),filters = 'hgnc_symbol', 
             values = gs,mart=ensembl)
eids=gsInfo[[2]]
names(eids)=gsInfo[[1]]
eids=eids[gs]
eids=eids[!is.na(eids)]
rma_nna=as.matrix(dat[names(eids),])
rownames(rma_nna) = eids
temp = gsva(rma_nna,canonicalC2BroadSets,min.sz=5,max.sz=500,rnaseq=T,mx.diff=TRUE,method="gsva",annotation="org.Hs.eg.db")
expr_gs = temp$es.obs
group <- factor(as.vector(meta[,"Drink"]))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_gs, design)
contrast.matrix <- makeContrasts("Y-N", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
allGeneSets <- topTable(fit2,number=nrow(fit2), adjust="BH",coef="Y-N")
PvalueCutoff <- 0.05
logFCcutoff <- log2(1.3)
DEgeneSets <- topTable(fit2, coef="Y-N", number=Inf, p.value=PvalueCutoff,adjust="none")
DEgeneSets
GSVAsco <- expr_gs[rownames(DEgeneSets), ]
pdf("GSVA-DEGeneSets-drinking vs. non-drinking.pdf", 18,12)
pheatmap(GSVAsco,cluster_rows=T,cluster_cols=T, treeheight_row=15,treeheight_col = 15, annotation_col= meta,
         scale = "row", annotation_legend=T,legend=T,show_rownames = T,show_colnames = T,
         cellheight=10, cellwidth=10,fontsize=8, fontsize_row = 8, 
         main = "Unadjusted P-value < 0.05",
         color = colorRampPalette(c("white", "yellow", "red"))(100))
dev.off()

#expression of genes correlated with survival (TCGA-KIRC)
gene_co = read.table("D:/Crownbio/Projects/Zhangning kidney/606 genes correlate with survival_TCGA.txt", h=T, sep = "\t",
                     stringsAsFactors = F, row.names = 1, check.names = F)
samples = row.names(subset(anno, Group == "Tumor" & Type == "Clear_cell"))
pdf("Heatmap-606 genes correlated with KIRC patient survival.pdf",16,12)
pheatmap(expr[row.names(gene_co),samples],cluster_rows=F,cluster_cols=T,treeheight_col = 10, scale = "row", 
         annotation_col = anno, annotation_row = gene_co, annotation_legend = T,legend=T,
         show_rownames = F,show_colnames = T, fontsize=6,  color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()




#check expression of mTOR signaling pathways
eids = canonicalC2BroadSets[[row.names(DEgeneSets)[1]]]@geneIds
gsInfo=getBM(attributes=c('hgnc_symbol','entrezgene','gene_biotype'),filters = 'entrezgene', 
             values = eids,mart=ensembl)
genes = gsInfo[[1]]
fit <- lmFit(log2(expr_g+1), design)
fit2 <- eBayes(fit)
results <- topTable(fit2,number=nrow(fit2),adjust.method = "BH",coef="groupTumor")
pheatmap(results[genes,1,drop=F],cluster_rows=F,cluster_cols=F, treeheight_row=15,treeheight_col = 15, 
         scale = "none", annotation_legend=T,legend=T,show_rownames = T,show_colnames = T,
         cellheight=15, cellwidth=15,fontsize=12, fontsize_row = 12, 
         color = colorRampPalette(c( "white", "red"))(100))



################################################################################
#EMT
markers = c("CDH1","CDH2","VIM","ACTA1","FN1","CTNNB1","SNAI1","TWIST1","SNAI2")
marker_expr = expr_g[markers,cc]
sample_anno = anno[colnames(expr_g),c(1:6,8,10,11)]
library("pheatmap")
pheatmap(marker_expr,cluster_rows=F,cluster_cols=T,treeheight_col = 15,scale = "row", annotation_col = sample_anno,
         annotation_legend = T,legend=T, show_rownames = T,show_colnames = T,
         fontsize=11,  color = colorRampPalette(c("blue", "white", "red"))(100))

marker_expr = as.data.frame(t(marker_expr))
marker_expr[,"Stage"] = as.factor(anno[rownames(marker_expr),"Stage"])
marker_expr[,"Grade"] = as.factor(anno[rownames(marker_expr),"Grade"])


proportion=table(marker_expr[,"Stage"])/nrow(marker_expr)
pdf("Boxplot of EMT expression grouped by stage.pdf",14,10)
par(mfrow = c(3, 3))
library(ggplot2)
library(ggsignif)
for(i in markers){
  test = kruskal.test(as.formula(paste(i, "Stage",sep = "~")), data = marker_expr)
  pv = test$p.value
  boxplot(as.formula(paste(i, "Stage",sep = "~")), data = marker_expr, width=proportion, at =c(1.3,2.3,3.3,4.3,5.3),
          range=0, whisklty = 0, staplelty = 0,xlab = "Stage", ylab = paste(i,"expression",sep = " "), 
          main = paste("p = ",pv,sep = ""))
  stripchart(as.formula(paste(i, "Stage",sep = "~")), data = marker_expr, vertical=TRUE, method="stack",add=TRUE,
             at =c(1.3,2.3,3.3,4.3,5.3), pch=19,col = c("red"))
}
dev.off() 

proportion=table(marker_expr[,"Grade"])/nrow(marker_expr)
pdf("Boxplot of EMT expression grouped by grade.pdf",14,10)
par(mfrow = c(3, 3))
for(i in markers){
  test = kruskal.test(as.formula(paste(i, "Grade",sep = "~")), data = marker_expr)
  pv = test$p.value
  boxplot(as.formula(paste(i, "Grade",sep = "~")), data = marker_expr, width=proportion, at =c(1.2,2.2,3.2,4.2),
          range=0, whisklty = 0, staplelty = 0,xlab = "Grade", ylab = paste(i,"expression",sep = " "), 
          main = paste("p = ",pv,sep = ""))
  stripchart(as.formula(paste(i, "Grade",sep = "~")), data = marker_expr, vertical=TRUE, method="stack",add=TRUE,
             at =c(1.2,2.2,3.2,4.2), pch=19,col = c("red"))
}
dev.off() 



#immunophenotyping
###############################################################################
markers = read.table("D:/Crownbio/Projects/Zhangning kidney/Immune_markers.txt", h=T, sep = "\t",
                     stringsAsFactors = F, check.names = F)
marker_expr = expr_g[markers[,1],]
sample_anno = anno[colnames(marker_expr),c(3:6,23,24,8:10)]
pdf("Immunophenotyping of ccRCC.pdf",20,15)
library("pheatmap")
pheatmap(marker_expr,cluster_rows=F,cluster_cols=T,treeheight_col = 15,scale = "row", annotation_col = sample_anno,
         annotation_legend = T,legend=T, show_rownames = T,show_colnames = T,
         fontsize=11,  color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()
group1 = c("ZN64T","ZN30T","ZN37T","ZN10","ZN18","ZN13","ZN21","ZN69","ZN62T","ZN29T","ZN14",
           "ZN2","ZN17","ZN4","ZN81","ZN80","ZN82")
group2 = c("ZN47T","ZN5","ZN22","ZN72","ZN41T","ZN73","ZN12","ZN66T","ZN74","ZN50T","ZN36","ZN31",
           "ZN63T","ZN54T","ZN53T","ZN78","ZN77","ZN58T","ZN68","ZN35T","ZN34T","ZN1","ZN65T","ZN46T",
           "ZN42T","ZN71","ZN60T","ZN57T","ZN16","ZN49T","ZN59T","ZN38T","ZN45T","ZN39T")
normal= c("ZN45N","ZN39N","ZN30N","ZN35N","ZN42N","ZN41N","ZN37N","ZN34N","ZN29N","ZN38N","ZN53N")
group3 = c("ZN26","ZN75","ZN56T","ZN61T")
ordered=c(group1,group2,normal,group3)
pdf("Immunophenotyping of ccRCC ordered.pdf",20,15)
pheatmap(marker_expr[,rev(ordered)],cluster_rows=F,cluster_cols=F,treeheight_col = 15,scale = "row", annotation_col = sample_anno,
         annotation_legend = T,legend=T, show_rownames = T,show_colnames = T,
         fontsize=11,  color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# immunophenotyping for TCGA tumor data
TCGA_T_g = TCGA_T
TCGA_T_g[,"Gene"] = sub("\\|.+", "", row.names(TCGA_T_g))
TCGA_T_g = TCGA_T_g[!duplicated(TCGA_T_g[,"Gene"]),]
row.names(TCGA_T_g) = TCGA_T_g[,"Gene"]
TCGA_T_g[,"Gene"] = NULL

marker_expr = TCGA_KIRC_g[markers[,1],]
sample_anno = classes[colnames(TCGA_KIRC_g),]
library("pheatmap")
pdf("Immunophenotyping of TCGA-KIRC.pdf",20,15)
hm <- pheatmap(marker_expr,cluster_rows=F,cluster_cols=T,treeheight_col = 15,scale = "row", annotation_col = sample_anno,
               annotation_legend = T,legend=T, show_rownames = T,show_colnames = T,
               fontsize=3,  color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

order <- colnames(marker_expr)[hm$tree_col$order] #obtain order of samples in heatmap hierarchy 
classes_immune=read.table("KIRC-immuno classification.txt",sep = "\t",h=T,stringsAsFactors = F,check.names = F,row.names = 1)

#Survival analysis in immuno-classification
library("survival")
library("coin")
temp = merge(clinical_KIRC, classes_immune, by.x= 0, by.y=0,all.y=T)
temp[,"merge"]=paste(temp[,"group"],temp[,"Immunophenotyping"],sep = "_")
#temp = na.omit(temp)
surv <- survfit(Surv(time, vital_status) ~ Immunophenotyping, type="kaplan-meier", conf.type="log", data = temp)
print(surv)
fit <- survdiff(Surv(time, vital_status) ~ Immunophenotyping, data = temp)
pv <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
pdf("Survival analysis on KIRC immunophenotyping.pdf",12,8)
plot(survfit(Surv(time, vital_status) ~ Immunophenotyping, type="kaplan-meier", conf.type="log",data = temp), main = paste("p =",as.character(pv),sep = " "),
     col=c("red","brown","blue"), ylab = "Probability", xlab = "Survival Time in Days",mark.time=TRUE, mark=3)
legend("topright", legend = c("active", "inactive","tolerant"),col=c("red","brown","blue"),bty = "n",title = "Immunophenotyping",lty=1)
dev.off()

library(survival) #Cox's proportional hazards model
fit <- coxph(Surv(time, vital_status) ~ Immunophenotyping + group + race + Stage_group + Grade_group + Metastatic, data = temp)
plot( survfit( fit),xlab="Time (in months)", ylab="Survival Function")
summary(fit)

surv <- survfit(Surv(time, vital_status) ~ merge,type="kaplan-meier", conf.type="log", data = temp)
print(surv)
fit <- survdiff(Surv(time, vital_status) ~ merge, data = temp)
pv <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
pdf("Survival analysis on KIRC immunophenotyping and classification.pdf",12,8)
plot(survfit(Surv(time, vital_status) ~ merge, type="kaplan-meier", conf.type="log",data = temp), main = paste("p =",as.character(pv),sep = " "),
     col=c("black","red","brown","yellow","green","darkgreen","darkblue","blue","Cyan"), ylab = "Probability", 
     xlab = "Survival Time in Days",mark.time=TRUE, mark=3)
legend("topright", legend = c("Class1_active", "Class1_inactive","Class1_tolerant","Class2_active", "Class2_inactive",
                              "Class2_tolerant","Class3_active", "Class3_inactive","Class3_tolerant"),
       col=c("black","red","brown","yellow","green","darkgreen","darkblue","blue","Cyan"),bty = "n",title = "Class",lty=1)
dev.off()
#######################################################################################

#EPIC
##################################################################################
library("EPIC")
score = EPIC(expr_g[,rev(ordered)], reference = NULL, mRNA_cell = NULL, mRNA_cell_sub = NULL,
     sigGenes = NULL, scaleExprs = TRUE, withOtherCells = TRUE,
     constrainedSum = TRUE, rangeBasedOptim = FALSE)
cells = c("Bcells", "CD4_Tcells", "CD8_Tcells","Macrophages", "NKcells","CAFs","Endothelial","otherCells")
fraction = t(score$cellFractions[,cells])
sample_anno = anno[colnames(fraction),c(3:6)]
write.table(fraction,"Cell fractions for ccRCC.txt",sep = "\t",row.names = T)

#estimate the correlation between TMB and TME
TMB=read.table("Number of somatic mutations for paired samples.txt",sep = "\t",header = T,row.names = 1)
a=merge(TMB,t(fraction),by.x=0,by.y=0)
row.names(a)=a[,1]
a[,1]=NULL
write.table(a,"TMB and cell fractions for 11 paired ccRCC.txt",sep = "\t",row.names = T)

b=data.frame(matrix(NA,nrow = ncol(a)-1, ncol = 2))
row.names(b)=colnames(a)[2:ncol(a)]
colnames(b)=c("Pearson correlation to TMB","P-value")
for(i in 2:ncol(a)){
  x=cor.test(a[,1],a[,i])
  b[i-1,]=c(cor(a[,1],a[,i]),x$p.value)
}

#stacked bar plots
pdf("Relative fraction of stromal cell types by EPIC.pdf ",12,8)
par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
colors=c("green","yellow","orange","red","purple","cyan","blue", "grey")
barplot(fraction, col=colors, width=1,space = 0,las=2, xlab = "Sample",ylab="Relative proportion",
        cex.axis = 0.6,cex.names = 0.6)
legend("topright", inset=c(-0.15,0),fill=colors, legend=row.names(fraction),cex = 0.6)
dev.off()

#vocano plot for immune active vs. other tumors
vocano = data.frame(matrix(NA, nrow = nrow(fraction), ncol = 2))
row.names(vocano) = row.names(fraction)
colnames(vocano) = c("log mean ratio","p-value")
library("gtools")
for(i in 1:nrow(fraction)){
  mean1 = mean(fraction[i,c(group1,group2)])
  mean2 = mean(fraction[i,c(normal,group3)])
  vocano[i,1] = log(mean1/mean2)
  vocano[i,2] = t.test(fraction[i,c(group1,group2)], fraction[i,c(normal,group3)], paired = F)$p.value
}
emf("vocano plot EPIC active tolerant vs. inactive normal.emf")
plot(vocano[,1],vocano[,2],log="y",xlab="log10 mean ratio",ylab="p-value",xlim=c(-10, 10),ylim=c(1,0.0000000001))
text(vocano[,1],vocano[,2],labels = row.names(vocano),pos=3,cex=0.9)
dev.off()


TCGA_KIRC_g=TCGA_KIRC
TCGA_KIRC_g[,"Gene"] = sub("\\|.+","",row.names(TCGA_KIRC_g))
TCGA_KIRC_g = TCGA_KIRC_g[!duplicated(TCGA_KIRC_g[,"Gene"]),]
row.names(TCGA_KIRC_g) = TCGA_KIRC_g[,"Gene"]
TCGA_KIRC_g[,"Gene"] = NULL
marker_expr = TCGA_KIRC_g[markers[,1],]
sample_anno = data.frame(matrix(NA,nrow = ncol(marker_expr),ncol=1))
row.names(sample_anno) = colnames(marker_expr)
pheatmap(marker_expr,cluster_rows=F,cluster_cols=T,treeheight_col = 30,scale = "row",
         annotation_legend = T,legend=T, show_rownames = T,show_colnames = F,
         fontsize=11,  color = colorRampPalette(c("blue", "white", "red"))(100))

################################################################################

#Estimate tumor purity
################################################################################
sc_expr = read.table("C:/Linda/Crownbio/Projects/Zhangning kidney/ccRCC/manuscript/Frontiers in Oncology/re-analysis/GSM1887306_Pt_mRCC_SC_01.TPM.txt",sep = "\t",header = T)
sc_expr = sc_expr[,c(1,7)]
temp = read.table("C:/Linda/Crownbio/Projects/Zhangning kidney/ccRCC/manuscript/Frontiers in Oncology/re-analysis/GSM1887312_Pt_mRCC_SC_10.TPM.txt",sep = "\t",header = T)
temp = temp[,c(1,7)]
sc_expr = merge(sc_expr,temp,by.x=1, by.y =1)
temp = read.table("C:/Linda/Crownbio/Projects/Zhangning kidney/ccRCC/manuscript/Frontiers in Oncology/re-analysis/GSM1887315_Pt_mRCC_SC_11.TPM.txt",sep = "\t",header = T)
temp = temp[,c(1,7)]
sc_expr = merge(sc_expr,temp,by.x=1, by.y =1)
temp = read.table("C:/Linda/Crownbio/Projects/Zhangning kidney/ccRCC/manuscript/Frontiers in Oncology/re-analysis/GSM1887316_Pt_mRCC_SC_12.TPM.txt",sep = "\t",header = T)
temp = temp[,c(1,7)]
sc_expr = merge(sc_expr,temp,by.x=1, by.y =1)
temp = read.table("C:/Linda/Crownbio/Projects/Zhangning kidney/ccRCC/manuscript/Frontiers in Oncology/re-analysis/GSM1887330_Pt_mRCC_SC_36.TPM.txt",sep = "\t",header = T)
temp = temp[,c(1,7)]
sc_expr = merge(sc_expr,temp,by.x=1, by.y =1)
row.names(sc_expr)=sc_expr[,1]
sc_expr[,1]=NULL
colnames(sc_expr) = c("SCtumor1","SCtumor2","SCtumor3","SCtumor4","SCtumor5")
row.names(sc_expr)=sub("\\.([0-9]{1,2})$","",row.names(sc_expr))

expr_fpkm = read.csv("C:/Linda/Crownbio/Projects/Zhangning kidney/ccRCC/manuscript/Frontiers in Oncology/re-analysis/RNAseq_expr_FPKM_11normals.csv",header = T)
ids = sub("^.*\\(","",expr_fpkm[,1])
ids = sub("\\)$","",ids)
row.names(expr_fpkm)=ids
symbols=sub("\\(.+\\)$","",expr_fpkm[,1])
expr_fpkm[,1]=symbols
expr_fpkm [,2:12]= 2^expr_fpkm[,2:12]

expr_merge=merge(expr_fpkm,sc_expr,by.x=0,by.y=0)
expr_merge = expr_merge[!duplicated(expr_merge[,2]),]
row.names(expr_merge)=expr_merge[,2]
expr_merge[,1:2] = NULL
write.table(expr_merge,"Reference expression matrix for ccRCC tumor and normal.txt",
            sep = "\t",row.names = T,quote = F)

library(ggplot2)
library("caret")
library("matrixStats")
genes <- apply(expr_merge , 1, function(x) (length(x[x>1])/length(x) >=0.2)) #in one row, at least 10% samples have expression level >0
expr_filtered = expr_merge[genes,]
nzv <- nearZeroVar(t(expr_filtered), saveMetrics= TRUE)
genes = row.names(nzv[nzv[,4]=="FALSE",])
expr_filtered = expr_filtered[genes,]

library("limma")
sample=colnames(expr_merge)
class=c(rep("Normal",11),rep("Tumor",5))
temp=cbind(class)
row.names(temp)=sample
group <- factor(as.vector(temp[,1]))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_filtered, design)
contrast.matrix <- makeContrasts("Normal-Tumor", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
all <- topTable(fit2,number=nrow(fit2), adjust="BH",coef="Normal-Tumor")
de=all[all$adj.P.Val<0.001,]
de = de[order(de[,1],decreasing = F),]
genes=row.names(de)[1:100]
de = de[order(de[,1],decreasing = T),]
genes=c(genes, row.names(de)[1:100])

x=cbind(rowMeans(expr_merge[,12:16]),rowMeans(expr_merge[,1:11]))
colnames(x)=c("Tumor cells","Normal cells")
ref = list()
ref[["refProfiles"]] = x
ref[["sigGenes"]] = genes

library("EPIC")
score = EPIC(expr_g, reference = ref, mRNA_cell = NULL, mRNA_cell_sub = NULL,
             sigGenes = NULL, scaleExprs = TRUE, withOtherCells = TRUE,
             constrainedSum = TRUE, rangeBasedOptim = FALSE)
score$mRNAProportions
tumor_purity=score$cellFractions
write.table(tumor_purity,"tumor purity estimation by EPIC.txt",sep = "\t",quote = F)

###########################################################################################3
#clustering CKC cc samples by NMF
dat = as.matrix(expr_t)
sample_anno = anno[colnames(dat),]
avg = mean(dat)
sd = sd(dat)
temp = (dat-avg)/sd
mad=rowMax(temp)-rowMeans(temp)
mad=mad[order(mad,decreasing = T)]

size=c(1500,2000,3000,4000,5000,6000,7000)
library("NMF")
for(k in size){
  top = names(mad[1:k])
  temp1 = temp[top,]
  temp1 = nneg(temp1)
  estim.r=nmf(temp1,rank=2:8,nrun=200,seed="random")
  plot(estim.r)
  pdf(paste("EstimateNMF-consensus matrices-cc-",as.character(k),"genes.pdf",sep = ""),20,13)
  consensusmap(estim.r, annCol = sample_anno[,c(1,3:6)],labCol = NA, labRow = NA)
  dev.off()
}

#rank=3
top = names(mad[1:3000])
temp1 = temp[top,]
temp1 = nneg(temp1)
res <- nmf(temp1, 3, nrun=200,seed="random")
cons <- consensusmap(res, annCol = sample_anno[,c(3:6,23,24)])
classes = lapply(cut(cons$Rowv,0.5)$lower, function(l) rapply(l, function(i) i))
class_anno = data.frame(matrix(NA,nrow = length(cc), ncol = 0))
row.names(class_anno)=cc
for(i in 1:length(classes)){
  for(j in 1:length(classes[[i]])){
    class_anno[classes[[i]][j],"NMFclass"]=paste("Class",as.character(i),sep = "")
  }
}

group <- factor(as.vector(class_anno[,1]))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(log2(dat[top,]+1), design)
contrast.matrix <- makeContrasts("Class1-Class2","Class1-Class3", "Class2-Class3",levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
de = data.frame()
for(k in colnames(contrast.matrix)){
  results =topTable(fit2,number=Inf, coef=k, p.value = 0.05,adjust="BH")
  results[,"Gene"] = row.names(results)
  de = rbind(de, results)
}
de = de[order(de[,5]),]
de = de[!duplicated(de[,7]),]
row.names(de) = 1:nrow(de)
de = de[,7]

results =topTable(fit2,number=Inf, coef="Class1-Class3", adjust="BH")
write.table(results, "DE Class1 vs Class3.txt",sep = "\t",quote = F)
results =topTable(fit2,number=Inf, coef="Class2-Class3", adjust="BH")
write.table(results, "DE Class2 vs Class3.txt",sep = "\t",quote = F)

sample_anno = merge(anno[colnames(dat),c(3:6,23,24)],class_anno,by.x=0,by.y =0)
row.names(sample_anno) = sample_anno[,1]
sample_anno[,1] = NULL
sample_anno = sample_anno[order(sample_anno[,7],sample_anno[,6]),]
ordered = row.names(sample_anno)
n = seq(41,100, by = 2)
accuracy = c()
for (i in n){
  temp1 = dat[de[1:i],ordered]
  temp2 = TCGA_T_g[de[1:i],]
  temp = cbind(temp1,temp2)
  cor=cor(temp,method="spearman")
  distance <- as.dist(sqrt(1 - cor^2))
  classes = as.data.frame(cutree(hclust(distance), k = 3)) 
  classes = merge(class_anno,classes, by.x = 0, by.y = 0)
  classes[,2] = sub("Class","",classes[,2])
  accuracy = c(accuracy,nrow(classes[classes[,2]==classes[,3],])/nrow(classes))
}
pdf("Accuracy for variety of DE genes.pdf")
plot(n, accuracy, xlabel="Number of DE genes", ylabel = "Accuracy")
dev.off()
#when n = 43, acurracy max = 0.9637
pdf("Heatmap-NMF classification cc with 43 DE genes.pdf",20,14)
library("pheatmap")
pheatmap(dat[de[1:43],ordered],cluster_rows=T,cluster_cols=F,treeheight_col = 15,scale = "row", 
         annotation_col = sample_anno, annotation_legend = T,legend=T, show_rownames = F,show_colnames = T,
         fontsize=11,  color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

temp1 = dat[de[1:43], ordered]
temp2 = TCGA_T_g[de[1:43],]
temp = cbind(temp1,temp2)
sample_anno = merge(pheno[colnames(temp),],class_anno,by.x=0,by.y =0, all.x=T)
row.names(sample_anno) = sample_anno[,1]
sample_anno[,1] = NULL
cor=cor(temp,method="spearman")
pdf("Heatmap-pairwise correlation of ccRCC.pdf",12,10)
pheatmap(cor,cluster_rows=T,cluster_cols=T,treeheight_col = 15,scale = "none", annotation_col = sample_anno,
         annotation_legend = T,legend=T, show_rownames = F,show_colnames = F,
         fontsize=11,  color = colorRampPalette(c("blue", "black", "yellow"))(100))
dev.off()
distance <- as.dist(sqrt(1 - cor^2))
plot(hclust(distance), main="Dissimilarity = 1 - Correlation", xlab="")
classes = as.data.frame(cutree(hclust(distance), k = 3)) 
classes = merge(classes,pheno[,c(2:7,9)],by.x=0,by.y=0)
classes[,2]= as.factor(paste("Class",classes[,2],sep = ""))
classes = classes[order(classes[,2],classes[,7],classes[,5],classes[,9]),]
write.table(classes,"ccRCC classification for CKC and TCGA KIRC-43 genes.txt",row.names = F,sep = "\t")
row.names(classes) = classes[,1]
classes[,1] = NULL
colnames(classes)[1] = "Class"

#ssGSVA on VEGF related genesets grouped by classification
##################################################################
library("GSVA")
temp= merge(dat[,ordered],TCGA_T_g,by.x=0,by.y=0)
row.names(temp) = temp[,1]
temp[,1] = NULL
gSet <- getGmt("D:/Crownbio/Projects/Zhangning kidney/VEGF related genesets.gmt")
results=gsva(as.matrix(temp),gSet,rnaseq=T,mx.diff=TRUE,method="ssgsea")
signature = row.names(results)
merged = merge(t(results),classes,by.x=0,by.y=0)
library("ggpubr")
for(i in 1:length(signature)){
  ggboxplot(merged, x = "Class", y = signature[i], color = "Class", palette = "jco",add = "jitter")
  + stat_compare_means()
  ggsave(paste("ssGSEA_for_classification_",signature[i],".png",sep = ""), plot = last_plot(),width = 12, height = 10)
}


#survival analysis grouped by classes
###################################################################
clinic_TCGA = anno_TCGA[colnames(TCGA_T),c(8:10,2:6,11:13)]
clinic_TCGA[] <- lapply(clinic_TCGA, as.character)
clinic_TCGA[clinic_TCGA[,1] == "alive",1] = 0
clinic_TCGA[clinic_TCGA[,1] == "dead",1] = 1
clinic_TCGA[,1:3] = apply(clinic_TCGA[,1:3], 2, function(x) as.numeric(x));
clinic_TCGA[,"time"] = clinic_TCGA[,2]
clinic_TCGA[is.na(clinic_TCGA[,2]),"time"]= clinic_TCGA[is.na(clinic_TCGA[,2]),3]
clinic_TCGA = clinic_TCGA[,c(1,12,4:11)]

samples=colnames(temp2)[colnames(temp2) %in% row.names(clinic_TCGA)]
classes = read.table("588 ccRCC classes.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,row.names = 1)
clinic_TCGA = merge(clinic_TCGA[samples,],classes[,1,drop=F],by.x=0, by.y=0)
row.names(clinic_TCGA) = clinic_TCGA[,1]
clinic_TCGA[,1]=NULL
clinic_TCGA[,11]= as.factor(clinic_TCGA[,11])
colnames(clinic_TCGA)[11] = "group"

TCGA_clustering=read.table("D:/Crownbio/Projects/Zhangning kidney/TCGA/TCGA-KIRC_2016cellreports_clustering_492.txt",
                                header = T,sep = "\t",stringsAsFactors = F,check.names = F,row.names = 1)
library("survival")
library("coin")
fit <- survdiff(Surv(time, vital_status) ~ group, data = clinic_TCGA)
pv <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
surv <- survfit(Surv(time, vital_status) ~ group, data = clinic_TCGA) 
print(surv)

temp = merge(clinic_TCGA, TCGA_clustering, by.x= 0, by.y=0)
#temp = na.omit(temp)
fit <- survdiff(Surv(time, vital_status) ~ group, data = temp)
pv <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
surv <- survfit(Surv(time, vital_status) ~ group, data = temp) 
print(surv)

surv <- survfit(Surv(time, vital_status) ~ Subtype, data = temp)
print(surv)
fit <- survdiff(Surv(time, vital_status) ~ Subtype, data = temp)
pv <- 1 - pchisq(fit$chisq, length(fit$n) - 1)



pdf("Survival analysis on KIRC classes.pdf",12,8)
plot(survfit(Surv(time, vital_status) ~ group, data = clinic_TCGA), main = paste("p =",as.character(pv),sep = " "),
     col=c("red","brown","blue"), ylab = "Probability", xlab = "Survival Time in Days",,mark.time=TRUE,mark=3)
legend("topright", legend = c("Class1", "Class2","Class3"),col=c("red","brown","blue"),bty = "n",title = "Class",lty=1)
dev.off()


###code ends#############

