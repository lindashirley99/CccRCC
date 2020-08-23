###This script is used to perform gene fusion analysis
###code starts#######################################

#summarize fusion data for all
path <- "D:/Crownbio/Projects/Zhangning kidney/fusion"
setwd(path)
pattern <- "N.filtered.fusion.list"
nfiles <- list.files(path = path, pattern = pattern, full.names = F)
normal = data.frame()
tumor = data.frame()
somatic = data.frame()
samples = c()
for(i in 1:length(nfiles)){
  sampname=sub("N.filtered.fusion.list","",nfiles[i])
  samples=c(samples, sampname)
  temp1 <- read.table(nfiles[i], h=T, sep = "\t", stringsAsFactors = F, check.names = F)
  for(j in 1:nrow(temp1)){
    x = c(temp1[j,2], temp1[j,6])
    sorted = sort(x)
    temp1[j,"Fus_Id"] = paste(sorted[1],sorted[2],sep = "_")
  }
  normal = rbind(normal, temp1)
  tfile = paste(sampname,"T.filtered.fusion.list",sep = "")
  temp2 <-  read.table(tfile, h=T, sep = "\t", stringsAsFactors = F, check.names = F)
  for(j in 1:nrow(temp2)){
    x = c(temp2[j,2], temp2[j,6])
    sorted = sort(x)
    temp2[j,"Fus_Id"] = paste(sorted[1],sorted[2],sep = "_")
  }
  tumor = rbind(tumor, temp2)
  temp = temp2[!temp2[,16] %in% temp1[,16],]
  somatic = rbind(somatic, temp)
}
write.table(somatic,file = paste(path,"12somatic.fusion.txt",sep = "/"), sep = "\t",row.names = F, col.names = T,quote = F)
library("plyr")
fus_freq=count(somatic[,16])
fus_freq=fus_freq[order(fus_freq[,2],decreasing = T),]
row.names(fus_freq) = 1:nrow(fus_freq)
pdf("Fusion summary-12 somatic.pdf",8,10)
barplot(fus_freq[1:14,2], width = 0.3, names.arg = fus_freq[1:14,1],main = "Fusion recurrence",col = rainbow(20),
        ylab = "Frequency in 12 samples",las=2,cex.axis = 0.5, cex.names = 0.5)
dev.off()

#all tumor samples
path = "/storage/wqian/clientdata/E3161-B1601zhangning/result/fusion/"
setwd(path)
pattern <- ".filtered.fusion.list"
files <- list.files(path = path, pattern = pattern, full.names = T)
pattern <- "N.filtered.fusion.list"
nfiles <- list.files(path = path, pattern = pattern, full.names = T)
files = files[!files %in% nfiles]
fusion = data.frame()
samples = c()
for(i in 1:length(files)){
  # sampname=sub(".filtered.fusion.list","",files[i])
  # samples=c(samples, sampname)
  temp <- read.table(files[i], h=T, sep = "\t", stringsAsFactors = F, check.names = F)
  temp = temp[!temp[,10]=="undetected",]
  temp = temp[!temp[,12]=="undetected",]
  #temp = temp[temp[,15]=="in-frame",]
  temp[,"Fus_Id"] = paste(temp[,2],temp[,6],sep = "_")
  temp=temp[!duplicated(temp[,16]),]
  fusion = rbind(fusion, temp)
}
links=data.frame(fusion[,c(3,5)],fusion[,5]+1,fusion[,c(7,9)],fusion[,9]+1,fusion[,16],fusion[,1])
colnames(links) = c("up_chr","up_start","up_end","dn_chr","dn_start","dn_end","fusion_ID","Sample")
write.table(links,file ="65tumors.fusion.filtered.txt", sep = "\t",row.names = F, col.names = T,quote = F)
write.table(fusion,file ="65tumors.fusion.filtered.txt", sep = "\t",row.names = F, 
            col.names = T,quote = F)
library("plyr")
fus_freq=count(fusion[,16])
fus_freq=fus_freq[order(fus_freq[,2],decreasing = T),]
row.names(fus_freq) = 1:nrow(fus_freq)
pdf("Fusion summary-65 tumors.pdf",8,10)
barplot(fus_freq[1:10,2], width = 0.3, names.arg = fus_freq[1:10,1],main = "Fusion recurrence",col = rainbow(20),
        ylab = "Frequency in 12 samples",las=2,cex.axis = 0.5, cex.names = 0.5)
dev.off()

fus_freq[,"Involved samples"] = ""
for(i in 1:nrow(fus_freq)){
  id <- fus_freq[i,1]
  for(j in 1:nrow(fusion)){
    if(id==fusion[j,16]){
      fus_freq[i,"Involved samples"] = paste(fus_freq[i,"Involved samples"], fusion[j,1], sep = ",")
    }
  }
}
fus_freq[,"Involved samples"] = sub("^,", "", fus_freq[,"Involved samples"])
temp = merge(fus_freq, fusion, by.x=1,by.y = 16,all.x = T)
temp=temp[!duplicated(temp[,1]),]
fus_freq = temp[,1:12]
write.table(fus_freq,file = paste(path,"65tumors fusion summary.txt",sep = "/"), sep = "\t", row.names = F, 
            col.names = T,quote = F)

fus_freq = fus_freq[order(fus_freq[,2],decreasing=T),]
for(i in 1:5){
  gene = fus_freq[i,"dw_gene"]
  samples = unlist(strsplit(fus_freq[i,"Involved samples"],split=","))
  temp = subset(anno,Group == "Tumor")
  sample_anno = data.frame(Sample = samples, fusion = "Y")
  sample_anno = rbind(sample_anno,data.frame(Sample = row.names(temp)[!row.names(temp) %in% samples], fusion = "N"))
  sample_anno = merge(sample_anno,log2(t(expr_g[gene,row.names(temp)])),by.x=1,by.y=0)
  library("ggpubr")
  ggboxplot(sample_anno, x = "fusion", y = gene, color = "fusion", palette = "jco",add = "jitter") + stat_compare_means()
  ggsave(paste("Boxplot of",gene,"expression grouped by fusion.png",sep=" "), plot = last_plot(),width = 12, height = 10)
}

###code ends##################
