###This R script perform mutational analysis
###code starts###############


#convert ANNOVAR annotation to MAF format for MutSigCV
library("maftools")
pattern <- ".finished.table.filtered.table"
files <- list.files(path = getwd(), pattern = pattern, full.names = F)
tumor = data.frame()
for(i in 1:length(files)){
  temp <- read.table(files[i], h=T, sep = "\t", stringsAsFactors = F, check.names = F)
  tumor = rbind(tumor, temp)
}
temp = tumor[tumor[,14] == "." | as.numeric(tumor[,14]) < 0.005,]
temp = temp[temp[,15] == "." | as.numeric(temp[,15]) < 0.005,]
temp = temp[temp[,7] == "exonic",]
proteincoding = c("frameshift deletion","frameshift insertion","nonsynonymous SNV","stopgain","stoploss",
                  "nonframeshift deletion","nonframeshift insertion")
temp = temp[temp[,10] %in% proteincoding,]
#temp = temp[temp[,50] == "PASS",]
cancercons = read.table("./mutation/cancerconsensusgenelist.table", h=F, sep = "\t", stringsAsFactors = F, check.names = F)
cancercons = cancercons[,1]
tumor_filtered = temp[temp[,46] %in% cancercons,]


write.table(tumor_filtered,file = "55tumors_filtered_consensus.txt",sep = "\t",row.names = F, col.names = T,quote = F)
tumor.maf <- annovarToMaf(annovar = "55tumors_filtered_consensus.txt", Center = NULL, refBuild = 'hg19',
                          tsbCol = 'Sample', table = 'ensGene',MAFobj = T)
temp = tumor_filtered[,c(1,46)]
temp = temp[!duplicated(temp[,c(1,2)]),]
library("plyr")
temp = count(temp[,2])
write.table(temp,"mutation counts for 55 tumors.txt",sep = "\t",row.names = F, col.names = T,quote = F)


#compare with TCGA
tumor.mutload = tcgaCompare(maf = tumor.maf, cohortName = 'CKC')
plotmafSummary(maf = tumor.maf, file = "55tumors MAF summary.pdf",width = 14, height = 10, rmOutlier = TRUE, showBarcodes=TRUE, textSize=1,
               addStat = NULL, dashboard = TRUE)
tumor.titv = titv(maf = tumor.maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = tumor.titv)

driver.maf <- annovarToMaf(annovar = "55tumors_driver_mutation.txt", Center = NULL, refBuild = 'hg19',
                           tsbCol = 'Sample', table = 'ensGene',MAFobj = T)
plotmafSummary(maf = driver.maf, file = "55 tumors driver MAF summary.pdf",width = 14, height = 10, rmOutlier = TRUE, 
               showBarcodes=TRUE, textSize=1, addStat = NULL, dashboard = TRUE)

pdf("Top 26 significantly mutated genes-driver.pdf",12,8)
oncoplot(maf = driver.maf, top = 26, fontSize = 12)
dev.off()

dnmt3a.lpop = lollipopPlot(maf = driver.maf, gene = 'PBRM1', AACol = 'AAChange', showMutationRate = TRUE, 
                           domainLabelSize = 3, defaultYaxis = FALSE)
tumor.titv = titv(maf = driver.maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = tumor.titv)
#Detecting cancer driver genes based on positional clustering.
somatic.sig = oncodrive(maf = driver.maf, AACol = 'AAChange', minMut = 3, pvalMethod = 'zscore')
head(somatic.sig)
plotOncodrive(res = somatic.sig, fdrCutOff = 0.1, useFraction = TRUE)

temp = read.table("55tumors_driver_mutation.txt",header = T, sep = "\t",stringsAsFactors = F)
#temp = temp[,c(1,46)]
samples = unique(temp[temp[,46] == "VHL",1]) #VHL mut samples
samples = row.names(sample_anno)[!row.names(sample_anno) %in% samples] #VHL wt samples
temp = temp[temp[,1] %in% samples,]
temp = temp[!duplicated(temp[,c(1,46)]),]
library("plyr")
freq = count(temp[,46])

pathway = read.table("Angiogenesis.txt",header = T, sep = "\t")
genes = as.vector(pathway[,1])
samples = unique(temp[temp[,46] %in% genes,1])


#TCGA mutation data
TCGA_mut = read.table("./TCGA/KImut.table",header = T, sep = "\t",stringsAsFactors = F)
TCGA_mut = TCGA_mut[TCGA_mut[,1] %in% samples_T,]
samples = unique(TCGA_mut[TCGA_mut[,46] == "VHL",1]) #VHL mut samples
#samples = unique(TCGA_mut[,1])[!unique(TCGA_mut[,1]) %in% samples] #VHL wt samples
sample_anno = anno_TCGA[samples,]
library("plyr")
freq = count(sample_anno[,6])

temp = TCGA_mut[TCGA_mut[,1] %in% samples,]
temp = temp[!duplicated(temp[,c(1,46)]),]

pathway = read.table("Angiogenesis.txt",header = T, sep = "\t")
genes = as.vector(pathway[,1])
samples = unique(temp[temp[,46] %in% genes,1])


samples = c("ZN21","ZN42T","ZN47T","ZN59T","ZN71")
temp = subset(anno,Group == "Tumor" & Type =="Clear_cell")
sample_anno = data.frame(Sample = samples, PBRM1mut = "Mutant")
sample_anno = rbind(sample_anno,data.frame(Sample = row.names(temp)[!row.names(temp) %in% samples], PBRM1mut = "Wild Type"))
sample_anno = merge(sample_anno,log2(t(expr_g["HIF1A",row.names(temp)])),by.x=1,by.y=0)
library("ggpubr")
ggboxplot(sample_anno,x = "PBRM1mut", y = "HIF1A", color = "PBRM1mut", palette = "jco",add = "jitter") + stat_compare_means()

temp = read.table("mutation frequencies in ccRCC.txt", h=T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)
fisher.test(as.matrix(temp), simulate.p.value=T)


temp = read.table("mutation frequencies in ccRCC.txt", h=T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)
fisher.test(as.matrix(temp), simulate.p.value=T)

temp = read.table("driver mutation CccRCC vs TCGA.txt", h=T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)
for(g in row.names(temp)){
  c1 = round(temp[g,1]*55/100)
  c2 = 55-c1
  c3 = round(temp[g,2]*451/100)
  c4 = 451-c3
  test = fisher.test(matrix(c(c1, c2, c3, c4), byrow = TRUE, 2, 2))
  temp[g,"p-value"] = test$p.value
}

write.table(temp, file = "driver mutation CccRCC vs TCGA_1.txt", sep = "\t", quote = F,row.names = T)

###code ends###################
