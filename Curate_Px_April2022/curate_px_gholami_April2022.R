#### Global Proteome Analysis of the NCI-60 Cell Line Panel ####
## Trying to replicate what Gholami et al. did in https://www.sciencedirect.com/science/article/pii/S221112471300380X?via%3Dihub#bib42 ##

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
setwd(wd)

lutp = read.table('prot_nci60_mapping.csv',header=TRUE,sep=",", as.is = TRUE, quote="");

#> unique(lutp$Dataset)
# [1] "kinome"                         "full proteome of 59 cell lines" "deep proteome of 9 cell lines"

idx = duplicated(lutp$Pnumber)
lutp = lutp[!idx,]

d1=read.table('proteinGroups_esc.txt',header=TRUE,sep="\t", as.is = TRUE, quote="");

names(d1) <- gsub("\\.","_", names(d1))
names(d1) <- gsub("__","_", names(d1))
names(d1) <- gsub("__","_", names(d1))
nm=names(d1)

idx = grep("LFQ", names(d1))

nm = names(d1[1,idx])
tmp = sapply(strsplit(nm,"\\_"), "[[",3)

length(intersect(tmp, unique(lutp$Pnumber)))

d2  =d1 [,c(1,idx)]

names(d2)

row.names(d2) = d2$Protein_IDs
d2$Protein_IDs = NULL

head(d2)
names(d2) <- gsub("\\LFQ_intensity_","", names(d2))
head(d2)

#### CURATION ###

## remove all proteins with the "REV_" tag
idx = grep("REV_",  rownames(d2))
d2 = d2[-idx,]

## remove all proteins with the "CON_" tag in it
idx = grep("CON_", rownames(d2))
d2 = d2[-idx,]

head(lutp)

lutp$cell2 = gsub("-", "",lutp$Cellline)
lutp$cell2 = gsub(" ", "",lutp$cell2)

lutp$cell2 = toupper(lutp$cell2)
lutp <- lutp[lutp$Dataset == "full proteome of 59 cell lines",]

head(lutp)

lutp[lutp$cell2 == "MCF7/ADRR",]$cell2 <- "MCF7ADR"
lutp$cell2[which(lutp$cell2 == "786O")] <- "7860"

`%notin%` <- function(x,y) !(x %in% y) 

#d3 = d2

#d2=d3

## nin = c("P001903", "P001564")

## SR is 2x in the data set -> P001388, P001750
sum(d2$P001388 == 0) # 4536 --> 514 # 4803 --> 525
sum(d2$P001750 == 0) # 5331 # 5597 --> 517 instead of 525

# Wir nehmen von beiden Samples dasjenige mit den wenigeren fehlenden Werten. Bei den SR musste man eine Auswahl treffen,
# weil unser integrativer Ansatz genau eine Probe verlangt.

## HOP-92 is 2x in the data set -> P003204, P003364
sum(d2$P003204 == 0) # 4552 --> 514  # 4792 --> 525
sum(d2$P003364 == 0) # 4693 # 4941 --> 525

nin =  c("P001750","P003364") # masking these data sets for later removal

library(tidyverse)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

## PCA with with both SR samples
i1750 <- which(colnames(d2) == "P001750")
i1388 <- which(colnames(d2) == "P001388")
t_d2 <- t(d2)
t_d2 <- t_d2[,apply(t_d2, 2, var, na.rm=TRUE) != 0]
d2.PCA <- prcomp(t_d2, center = F, scale. = T)
groups <- rep("BASIC", 61)
groups[i1750] <- "P001750"
groups[i1388] <- "P001388"
pdf("PCA_Px_BASIC_SR.pdf")
ggbiplot(d2.PCA, labels = rownames(t_d2), var.axes = F, obs.scale = 1, vavar.scale = 1, groups = groups)
dev.off()

## PCA with with both HOP92 samples
i3364 <- which(colnames(d2) == "P003364")
i3204 <- which(colnames(d2) == "P003204")
t_d2 <- t(d2)
t_d2 <- t_d2[,apply(t_d2, 2, var, na.rm=TRUE) != 0]
d2.PCA <- prcomp(t_d2, center = F, scale. = T)
groups <- rep("BASIC", 61)
groups[i3364] <- "P003364"
groups[i3204] <- "P003204"
pdf("PCA_Px_BASIC_HOP92.pdf")
ggbiplot(d2.PCA, labels = rownames(t_d2), var.axes = F, obs.scale = 1, vavar.scale = 1, groups = groups)
dev.off()

xt = read.table('transcr_exp_design.tsv',header=TRUE,sep="\t", as.is = TRUE, quote="", skip = 0);
#xt$cancer =  sapply(strsplit(xt$cell_linet,"\\_"), "[[",1)

xt$clt =  sapply(strsplit(xt$cell_linet,"\\_"), "[[",2)

for ( i in 1:(dim(d2)[2])) # add abbreviations CO, BR...etc
{
  print(i)
  if(names(d2)[i] %notin% nin)
  {
    idx = which(lutp$Pnumber == names(d2)[i])
    
    idx = which(xt$clt == lutp$cell2[idx[1]]) # cancer type from array data
    if(idx>0)
    {
      names(d2)[i] = paste(xt$cancer2[idx] , xt$clt[idx], sep=".")
    }
  }#nin
}

d2[-which(names(d2) == "P001750")] -> d2
d2[-which(names(d2) == "P003364")] -> d2

names(d2)[which(names(d2) == "BR.MCF7ADR")] <- "OV.ADRRES"
#https://www.ncbi.nlm.nih.gov/pubmed/16504380
#https://www.hindawi.com/journals/bmri/2019/6146972/

names(d2)[which(names(d2) == "BR.MDAMB435")] <- "ME.MDAMB435"
#https://www.nature.com/articles/npjbcancer20152
#https://www.ncbi.nlm.nih.gov/pubmed/28940260

rownames(d2) <- gsub("__", "", rownames(d2))
## split by ";" and only take the first entry
rownames(d2) <- gsub(":", "", rownames(d2))

rownames(d2) <- as.character(rownames(d2))

head(rownames(d2))
get_first_name <- function(x) {
  return(strsplit(x, ";")[[1]][1])
}
new_rownames <- c()
for (i in 1:length(rownames(d2))) {
  new_rownames <- c(new_rownames, get_first_name(rownames(d2)[i]))
}
rownames(d2) <- new_rownames
head(rownames(d2))

write.table(x = d2, "proteinGroups_esc_polished.txt", col.names = T, row.names = T, sep = "\t")
d_log10 <- log10(d2)
write.table(x = d_log10, "proteinGroups_esc_polished_log10.txt", col.names = T, row.names = T, sep = "\t")


cancer.types <- unique(xt$cancer2)
df <- setNames(data.frame(matrix(ncol = 9, nrow = length(rownames(d2)))), cancer.types)

#head(df[rowSums(df) == 59,])
#df[rowSums(df) == 59,]

#length(grep("BR.", colnames(d2))) #5
#length(grep("CNS.", colnames(d2))) #6
#length(grep("CO.", colnames(d2))) #7
#length(grep("LEUK.", colnames(d2))) #6
#length(grep("ME.", colnames(d2))) #9
#length(grep("LC.", colnames(d2))) #9
#length(grep("OVAR.", colnames(d2))) #7
#length(grep("PROSTATE.", colnames(d2))) #2
#length(grep("RE\\.", colnames(d2))) #8

for (cancer.type in cancer.types) {
  if (cancer.type == "RE") {
    cancer.type <- "RE\\."
  }
  subs <- d2[,grep(cancer.type, colnames(d2))]
  if (cancer.type == "RE\\.") {
    cancer.type <- "RE"
  }
  df[, cancer.type] <- rowSums(subs > 0)
}

marked_lines <- c()
for (line in 1:dim(df)[[1]]) {
  l <- df[line,]
  if (rowSums(l) == 59) {
    marked_lines <- c(marked_lines, line)
  }
}

d_final <- d2[marked_lines,]

write.table(x = d_final, "proteinGroups_esc_polished_filtered.txt", col.names = T, row.names = T, sep = "\t")
d_final_log10 <- log10(d_final)
write.table(x = d_final_log10, "proteinGroups_esc_polished_filtered_log10.txt", col.names = T, row.names = T, sep = "\t")

# http://52.71.54.154/help/course-materials/2015/UseBioconductorFeb2015/A01.5_Annotation.html
# https://support.bioconductor.org/p/9136239/
# options(connectionObserver = NULL)
# library(org.Hs.eg.db)
# cols <- c("SYMBOL")
# d_log10_ipi_to_gene_symbol <- select(org.Hs.eg.db, keys=rownames(d_log10), columns = cols, keytype = "IPI")
# duplicated_ipi <- duplicated(d_log10_ipi_to_gene_symbol$IPI)
# d_log10_ipi_to_gene_symbol_uniq <- d_log10_ipi_to_gene_symbol[!duplicated_ipi,]
# stopifnot(all(d_log10_ipi_to_gene_symbol_uniq$IPI == rownames(d_log10)))
# new_rownames <- paste(d_log10_ipi_to_gene_symbol_uniq$IPI, d_log10_ipi_to_gene_symbol_uniq$SYMBOL, sep = "_")
# rownames(d_log10) <- new_rownames
# write.table(x = d_log10, "proteinGroups_esc_polished_log10_symboles.txt", col.names = T, row.names = T, sep = "\t")
# 
# d_final_log10_ipi_to_gene_symbol <- select(org.Hs.eg.db, keys=rownames(d_final_log10), columns = cols, keytype = "IPI")
# duplicated_ipi <- duplicated(d_final_log10_ipi_to_gene_symbol$IPI)
# d_final_log10_ipi_to_gene_symbol_uniq <- d_final_log10_ipi_to_gene_symbol[!duplicated_ipi,]
# stopifnot(all(d_final_log10_ipi_to_gene_symbol_uniq$IPI == rownames(d_final_log10)))
# new_rownames <- paste(d_final_log10_ipi_to_gene_symbol_uniq$IPI, d_final_log10_ipi_to_gene_symbol_uniq$SYMBOL, sep = "_")
# rownames(d_final_log10) <- new_rownames
# write.table(x = d_final_log10, "proteinGroups_esc_polished_filtered_log10_symboles.txt", col.names = T, row.names = T, sep = "\t")

#### ENDING CURATION ####

