#clean up
rm(list = ls(all = TRUE)) # clear all variables
graphics.off()

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

load("../recount2_tcga_data/rse_gene_bone_marrow.Rdata")
# load("../recount2_tcga_data/rse_gene_lymph_nodes.Rdata")
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

m1 = m[,grepl("gdc_cases.submitter_id|sample_type|primary_site|gender|race|disease_type|tumor|risk_factors",names(m))]
head(m1)
rm(m,new_name)

# remove NA values in Metadata by removing the actual samples
if (length(which(is.na(m1$gdc_cases.demographic.gender))) > 0) {
  rse_gene <- rse_gene[,-which(is.na(m1$gdc_cases.demographic.gender))]
  m1 <- m1[-which(is.na(m1$gdc_cases.demographic.gender)),]
}
rse_gene <- scale_counts(rse_gene)
# only keep "Primary Tumor" and "Solid Tissue Normal"
rse_gene <- rse_gene[, m1$gdc_cases.samples.sample_type == "Primary Tumor" | m1$gdc_cases.samples.sample_type == "Solid Tissue Normal"]
m1 <- m1[m1$gdc_cases.samples.sample_type == "Primary Tumor" | m1$gdc_cases.samples.sample_type == "Solid Tissue Normal",]
## Primary Tumor
sum(m1$gdc_cases.samples.sample_type == "Primary Tumor") # 48
## Solid Tissue Normal
sum(m1$gdc_cases.samples.sample_type == "Solid Tissue Normal") # 0
# As we have 0 normal tissues samples we have to
#### ABORT MISSION ####