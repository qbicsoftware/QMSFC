rm(list = ls(all = TRUE)) # clear all  variables
graphics.off()

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# install.packages("rstudioapi") # run this if it's your first time using it to install
library(rstudioapi) # load it
# the following line is for getting the path of your current open file
current_path <- getActiveDocumentContext()$path 
# The next line set the working directory to the relevant one:
setwd(dirname(current_path ))
# you can make sure you are in the right directory
print(getwd())
wd <- getwd()

library(recount)
library(reshape2)
library(ggplot2)
library(DESeq2)

load("../recount2_tcga_data_2021-22-04/rse_gene_colorectal.Rdata")
dim(colData(rse_gene))
head(assay(rse_gene))

head(colData(rse_gene)$gdc_cases.submitter_id)
rse_gene@colData@rownames <- paste(colData(rse_gene)$gdc_cases.submitter_id,"id",1:dim(colData(rse_gene))[1],sep="_")
head(rse_gene@colData@rownames)
new_name = data.frame(elementMetadata(rse_gene))
row.names(rse_gene) <- paste(new_name[,1],new_name[,3],sep="_")
head(assay(rse_gene))

#check coldata with useful info
m = data.frame(colData(rse_gene))
m = m[colSums(!is.na(m)) > 0]

m1 = m[,grepl("gdc_cases.submitter_id|sample_type|primary_site|gender|race|disease_type|tumor|risk_factors|tissue",names(m))]
head(m1)
rm(m,new_name)
rse_gene <- scale_counts(rse_gene)

# save.image("recount2_co_DESeq.Rdata")
  
# load(file = ".RData")
# remove NA values in Metadata by removing the actual samples
rse_gene <- rse_gene[,-which(is.na(m1$gdc_cases.demographic.gender))]
m1 <- m1[-which(is.na(m1$gdc_cases.demographic.gender)),]
# only keep "Primary Tumor" and "Solid Tissue Normal"
rse_gene <- rse_gene[, m1$gdc_cases.samples.sample_type == "Primary Tumor" | m1$gdc_cases.samples.sample_type == "Solid Tissue Normal"]
m1 <- m1[m1$gdc_cases.samples.sample_type == "Primary Tumor" | m1$gdc_cases.samples.sample_type == "Solid Tissue Normal",]
## Primary Tumor

## Solid Tissue Normal

# xml_tumor_tissue_site "Colon" only
colon_only <- m1$xml_tumor_tissue_site == "Colon"
colon_only[is.na(colon_only)] <- FALSE
rse_gene <- rse_gene[, colon_only]
m1 <- m1[colon_only,]
## Primary Tumor
sum(m1$gdc_cases.samples.sample_type == "Primary Tumor") # 500
## Solid Tissue Normal
sum(m1$gdc_cases.samples.sample_type == "Solid Tissue Normal") # 41

coSE <- DESeqDataSet(rse_gene, design = ~gdc_cases.demographic.gender + gdc_cases.demographic.race + gdc_cases.samples.sample_type)
# coSE <- DESeqDataSet(rse_gene, design = ~gdc_cases.samples.sample_type + gdc_cases.demographic.gender)
coSE$gdc_cases.samples.sample_type <- factor(coSE$gdc_cases.samples.sample_type, levels = c("Solid Tissue Normal","Primary Tumor")) 
coSE <- DESeq(coSE, parallel = F)

# save.image("recount2_co_DESeq.Rdata")

res05 <- results(coSE, alpha=0.05, lfcThreshold = 2, contrast = c("gdc_cases.samples.sample_type", "Primary Tumor", "Solid Tissue Normal"))
res05 <- subset(res05, select = -c(stat, lfcSE, baseMean))
res05 <- res05[!is.na(res05$padj),]
res05 <- res05[res05$padj < 0.05,] # 1900
gene_names <- c()
for (i in 1:length(res05$padj)) {
  ens <- rownames(res05)[i]
  gene_name <- strsplit(ens, "_")[[1]][2]
  gene_names <- c(gene_names, gene_name)
}
res05$gene_name <- gene_names
write.table(res05, "colon_sig_0.05_FC2.tsv", append = F, quote = F, sep = "\t", dec = ".")

res <- results(coSE, contrast = c("gdc_cases.samples.sample_type", "Primary Tumor", "Solid Tissue Normal"))
res <- subset(res, select = -c(stat, lfcSE, baseMean))
gene_names <- c()
for (i in 1:length(res$padj)) {
  ens <- rownames(res)[i]
  gene_name <- strsplit(ens, "_")[[1]][2]
  gene_names <- c(gene_names, gene_name)
}
res$gene_name <- gene_names
write.table(res, "colon_DESeq2.tsv", append = F, quote = F, sep = "\t", dec = ".")

vsd <- vst(coSE, blind=FALSE)
norm_vals <- assay(vsd)
write.table(norm_vals, "vst_transformed_colon.read.counts.tsv", append = F, quote = F, sep = "\t", dec = ".")

# save.image("recount2_co_done.Rdata")
  