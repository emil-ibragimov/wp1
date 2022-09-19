library(EnhancedVolcano)
library('dplyr')
library('DESeq2')
library('ggplot2')

#upload data
setwd('C:/Users/mdb241/Documents/RNAseq/')
gene_annot <- read.csv('pig_annotation.csv')
pheno <- read.table('wp1_individ.tsv', header = TRUE)
count_table <- read.csv('counts.csv')
# count_table
colnames(count_table) <- c('gene', substring(colnames(count_table), 59, 76)[-1])
rownames(count_table) <- count_table$gene
count_table <- count_table[,-1]
barcode_map = read.csv('barcode_map.csv', sep = ';', dec = ',')
barcode_map['full_barcode'] <- paste(barcode_map$lane_id, barcode_map$barcode_id, sep = '_')
barcode_map <- barcode_map %>%
  rename('ID' = 'Pig_ID')
coldata <- pheno[,c('Pig_ID', 'Group', 'Sex', 'Body_weight_start_kg', 'Feed_units_per_day', 'Daily_body_weight_gain_g')]
barcode_map$Pig_ID <- paste('sam', barcode_map$Pig_ID, sep = '')
pca <- read.csv('pca.csv')

colnames(count_table) <- barcode_map$Pig_ID
coldata$Pig_ID <- paste('sam', coldata$Pig_ID, sep = '')

coldata <- coldata[coldata$Pig_ID %in% colnames(count_table),]
coldata <- coldata[match(colnames(count_table), coldata$Pig_ID),]

#exclude genes which are zero in more than 20 samples
count_zeros <- function(a){
  return(sum(a == 0))
}
coldata$Sex <- as.factor(coldata$Sex)
coldata$Group <- as.factor(coldata$Group)
coldata$Daily_body_weight_gain_kg <- coldata$Daily_body_weight_gain_g / 1000
coldata$FCR <- coldata$Daily_body_weight_gain_kg/coldata$Feed_units_per_day
coldata <- coldata[!coldata$Pig_ID %in% c('sam2503', 'sam2842', 'sam2873', 'sam9672'),]
count_table <- count_table[,!colnames(count_table) %in% c('sam2503', 'sam2842', 'sam2873', 'sam9672')]
n_obs_all <- dim(count_table)[2]
simple_rfi <- read.csv('simple_rfi.csv')

coldata <- merge(coldata,simple_rfi)




count_table_withoutzeros <- count_table[apply(count_table, 1, count_zeros) < 0.2*n_obs_all,]
count_table_withoutzeros <- count_table_withoutzeros[rowSums(count_table_withoutzeros) > 0.8*n_obs_all,]

coldata <- coldata[match(colnames(count_table_withoutzeros), coldata$Pig_ID),]

coldata_m <- coldata[coldata$Sex == 1,]
coldata_f <- coldata[coldata$Sex == 0,]

count_table_m <- count_table[,colnames(count_table) %in% coldata_m$Pig_ID]
count_table_f <- count_table[,colnames(count_table) %in% coldata_f$Pig_ID]



n_obs_m <- dim(count_table_m)[2]
count_table_withoutzeros_m <- count_table_m[apply(count_table_m, 1, count_zeros) < 0.2*n_obs_m,]
count_table_withoutzeros_m <- count_table_withoutzeros_m[rowSums(count_table_withoutzeros_m) > 0.8*n_obs_m,]


n_obs_f <- dim(count_table_f)[2]
count_table_withoutzeros_f <- count_table_f[apply(count_table_f, 1, count_zeros) < 0.2*n_obs_f,]
count_table_withoutzeros_f <- count_table_withoutzeros_f[rowSums(count_table_withoutzeros_f) > 0.8*n_obs_f,]

coldata$cc_rfi <- as.factor(ifelse(coldata$simple_RFI >= 0, 1, 0))
coldata_m$cc_rfi <- as.factor(ifelse(coldata_m$simple_RFI >= 0, 1, 0))
coldata_f$cc_rfi <- as.factor(ifelse(coldata_f$simple_RFI >= 0, 1, 0))

colnames(count_table_withoutzeros_m)
coldata_m$Pig_ID
colnames(coldata)
coldata$Body_weight_start_kg


dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + Body_weight_start_kg + FCR)

dds <- DESeq(dds, parallel = TRUE)

res_fcr <- results(dds)
res_sw <- results(dds, name = 'Body_weight_start_kg')
hist(res_sw$pvalue, breaks = 100)



dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + Body_weight_start_kg + Daily_body_weight_gain_kg + Feed_units_per_day)


dds <- DESeq(dds, parallel = TRUE)
hist(results(dds)$pvalue, breaks = 100)
res_full <- results(dds)
sum(res_full$padj<=0.1, na.rm = TRUE)



dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + sv1 + sv2 + Body_weight_start_kg + Daily_body_weight_gain_kg + Feed_units_per_day)


#full model
dds <- DESeq(dds, parallel = TRUE)
hist(results(dds)$pvalue, breaks = 100)
res_full_sv <- results(dds)
sum(res_full$padj<=0.1, na.rm = TRUE)
sum(res_full_sv$padj<=0.1, na.rm = TRUE)
write.csv(res_full_sv, 'res_full_sv.csv')


#start weight
mod = model.matrix(~Sex + Body_weight_start_kg, data=pheno)
mod0 = model.matrix(~Sex, data=pheno)
svobj = sva(edata,mod,mod0,n.sv=2)
sv_df <- data.frame(svobj$sv)
colnames(sv_df) <- c('sv1', 'sv2')
coldata <- coldata[,1:10]
coldata <- cbind(coldata, sv_df)

dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + sv1 + sv2 + Body_weight_start_kg)


dds <- DESeq(dds, parallel = TRUE)

res_sw_sv <- results(dds)
write.csv(res_sw_sv, 'res_sw_sv.csv')

#daily gain
mod = model.matrix(~Sex + Daily_body_weight_gain_kg, data=pheno)
mod0 = model.matrix(~Sex, data=pheno)
svobj = sva(edata,mod,mod0,n.sv=2)
sv_df <- data.frame(svobj$sv)
colnames(sv_df) <- c('sv1', 'sv2')
coldata <- coldata[,1:10]
coldata <- cbind(coldata, sv_df)

dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + sv1 + sv2 + Daily_body_weight_gain_kg)


dds <- DESeq(dds, parallel = TRUE)

res_dg_sv <- results(dds)
write.csv(res_dg_sv, 'res_dg_sv.csv')

#feed intake
mod = model.matrix(~Sex + Feed_units_per_day, data=pheno)
mod0 = model.matrix(~Sex, data=pheno)
svobj = sva(edata,mod,mod0,n.sv=2)
sv_df <- data.frame(svobj$sv)
colnames(sv_df) <- c('sv1', 'sv2')
coldata <- coldata[,1:10]
coldata <- cbind(coldata, sv_df)

dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + sv1 + sv2 + Feed_units_per_day)



dds <- DESeq(dds, parallel = TRUE)

res_fi_sv <- results(dds)
write.csv(res_fi_sv, 'res_fi_sv.csv')
to_plot <- plotCounts(dds, gene = 'ENSSSCG00000032434', intgroup = 'Feed_units_per_day', returnData = TRUE)
ggplot(to_plot, aes(Feed_units_per_day, count))+
  geom_point()
to_plot <- plotCounts(dds, gene = 'ENSSSCG00000032434', intgroup = 'FCR', returnData = TRUE)
ggplot(to_plot, aes(FCR, count))+
  geom_point()



#INCLUDE START WEIGHT!!!!!!!!!!!

#daily gain
mod = model.matrix(~Sex + Body_weight_start_kg + Daily_body_weight_gain_kg, data=pheno)
mod0 = model.matrix(~Sex + Body_weight_start_kg, data=pheno)
svobj = sva(edata,mod,mod0,n.sv=2)
sv_df <- data.frame(svobj$sv)
colnames(sv_df) <- c('sv1', 'sv2')
coldata <- coldata[,1:10]
coldata <- cbind(coldata, sv_df)

dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + sv1 + sv2 +Body_weight_start_kg + Daily_body_weight_gain_kg)


dds <- DESeq(dds, parallel = TRUE)

res_dg_sv <- results(dds)
write.csv(res_dg_sv, 'res_dg_sv_sw.csv')

#feed intake
mod = model.matrix(~Sex + Body_weight_start_kg + Feed_units_per_day, data=pheno)
mod0 = model.matrix(~Sex + Body_weight_start_kg, data=pheno)
svobj = sva(edata,mod,mod0,n.sv=2)
sv_df <- data.frame(svobj$sv)
colnames(sv_df) <- c('sv1', 'sv2')
coldata <- coldata[,1:10]
coldata <- cbind(coldata, sv_df)

dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + sv1 + sv2 + Body_weight_start_kg + Feed_units_per_day)



dds <- DESeq(dds, parallel = TRUE)

res_fi_sv <- results(dds)
write.csv(res_fi_sv, 'res_fi_sv_sw.csv')
#PERMUTED FEED INTAKE
library(gtools)
coldata$perm_Feed_units_per_day <- permute(coldata$Feed_units_per_day)
pheno$perm_Feed_units_per_day <- coldata$perm_Feed_units_per_day
mod = model.matrix(~Sex + perm_Feed_units_per_day, data=pheno)
mod0 = model.matrix(~Sex, data=pheno)
svobj = sva(edata,mod,mod0,n.sv=2)
sv_df <- data.frame(svobj$sv)
colnames(sv_df) <- c('sv1', 'sv2')
coldata <- coldata[,1:10]
coldata <- cbind(coldata, sv_df)
coldata$perm_Feed_units_per_day <- pheno$perm_Feed_units_per_day
dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + sv1 + sv2 + perm_Feed_units_per_day)



dds <- DESeq(dds, parallel = TRUE)

res_fi_per_sv <- results(dds)
sum(res_fi_sv$padj < 0.1, na.rm = TRUE)
hist(res_fi_per_sv$pvalue, breaks = 100)
write.csv(res_fi_sv, 'res_fi_sv.csv')



#CUT
coldata$cut4_Feed_units_per_day <- cut(coldata$Feed_units_per_day, breaks = 4)
coldata$cut4_Body_weight_start_kg <- cut(coldata$Body_weight_start_kg, breaks = 4)
coldata$cut4_Daily_body_weight_gain_kg <- cut(coldata$Daily_body_weight_gain_kg, breaks = 4)
coldata$cut4_rfi <- cut(coldata$simple_RFI, breaks = 4)
coldata$cut4_fcr <- cut(coldata$FCR, breaks = 4)

plot(coldata$cut4_Body_weight_start_kg, coldata$Body_weight_start_kg)
plot(coldata$cut4_Feed_units_per_day, coldata$Feed_units_per_day)

dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + cut4_Body_weight_start_kg + Feed_units_per_day)



dds <- DESeq(dds, parallel = TRUE)
results(dds)
resultsNames(dds)
res_fi_cut_sw <- results(dds)
sum(res_fi_cut_sw$padj < 0.1, na.rm = TRUE)
hist(res_fi_cut_sw$pvalue, breaks = 100)
#




#case control for feed intake
coldata$cc_fi <- as.factor(ifelse(coldata$Feed_units_per_day >= median(coldata$Feed_units_per_day), 1, 0))

dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + cc_fi)

dds <- DESeq(dds, parallel = TRUE)
hist(results(dds)$pvalue, breaks = 100)
res_cc_fi <- results(dds)

write.csv(res_cc_fi, 'res_cc_fi.csv')

#



#



columns(coldata$Body_weight_start_kg)
dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + Body_weight_start_kg)


hist(res_fcr$pvalue, breaks = 100)
resultsNames(dds)







dds <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros,
                              colData = coldata,
                              design= ~ Sex + Body_weight_start_kg + FCR)
dds <- DESeq(dds, parallel = TRUE)
hist(results(dds)$pvalue, breaks = 100)
res_fcr <- results(dds)






dds_m <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros_m,
                              colData = coldata_m,
                              design= ~ cc_rfi)

dds_m <- DESeq(dds_m, parallel = TRUE)


dds_f <- DESeqDataSetFromMatrix(countData = count_table_withoutzeros_f,
                                colData = coldata_f,
                                design= ~ cc_rfi)

dds_f <- DESeq(dds_f, parallel = TRUE)


table(coldata_m$cc_rfi)
table(coldata_f$cc_rfi)

dim(count_table_withoutzeros)


hist(results(dds)$pvalue, breaks = 100)
hist(results(dds_m)$pvalue, breaks = 100)
hist(results(dds_f)$pvalue, breaks = 100)
write.csv(results(dds), 'res_cor_all.csv')
write.csv(results(dds_m), 'res_cor_m.csv')
write.csv(results(dds_f), 'res_cor_f.csv')

design(dds)<-~Sex+Sex:cc_rfi
dds <- DESeq(dds, parallel = TRUE)
res_tmp1 <- results(dds, contrast=list('Sex0.cc_rfi1'))
res_tmp2 <- results(dds, contrast=list('Sex1.cc_rfi1'))
design(dds)<-~Sex+cc_rfi+Sex:cc_rfi
dds <- DESeq(dds, parallel = TRUE)
resultsNames(dds)
write.csv(results(dds,contrast = list('cc_rfi_1_vs_0')), 'res_cor_f2.csv')
write.csv(results(dds,contrast = list(c('cc_rfi_1_vs_0', 'Sex1.cc_rfi1'))), 'res_cor_m2.csv')

vsd <- vst(dds, blind = TRUE)
write.csv(assay(vsd), 'vst_all.csv')

results(dds)



vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="Sex")
dds
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IHW")

library("IHW")



resIHW <- results(dds, filterFun=ihw)
res <- results(dds)
sum(resIHW$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.1, na.rm=TRUE)
head(res[order(res$pvalue),],10)








head(coldata)
design(dds) <- ~Sex+cc_rfi
dds <- DESeq(dds, parallel = TRUE)
res_dg <- results(dds)

hist(res_dg$pvalue, breaks = 100)
hist(res_fi$pvalue, breaks = 100)
design(dds) <- ~Sex+Feed_units_per_day
dds <- DESeq(dds, parallel = TRUE)
res_fi <- results(dds)
head(coldata$cc_rfi)


res_fi <- read.csv('res_fi_sv.csv', row.names = 1)
res_fi <- res_fi[!is.na(res_fi$padj),]
de_genes<-rownames(res_fi)[res_fi$padj < 0.1]
write.csv(assay(dds)[de_genes,], 'norm_de.csv')

