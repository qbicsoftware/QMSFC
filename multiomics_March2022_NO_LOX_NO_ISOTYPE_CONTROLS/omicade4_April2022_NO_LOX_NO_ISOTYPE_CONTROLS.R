### R code from vignette source 'omicade4.Rnw'

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
#### --------------FACS ####

setwd(paste (wd,'/lists_tsv/',sep=""))    # assumes this subfolder exists and contains the protein summaries

files_m <- dir()    # get the filenames

BG = NULL
for (i in 1:length(files_m))   # loop through files
{
  fn = files_m[i];
  print(fn)
  d1 = read.table(fn,  header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,
                  stringsAsFactors=F,dec=".",fill=TRUE, as.is = T, check.names=FALSE)
  d1$srt = d1$WID = d1$receptor = NULL
  
  if(i==1)
  {BG = d1}
  else{
    BG = merge(BG,d1, x.by = "receptor_esc", y.by ="receptor_esc")
  }
  
  
}

setwd(wd)

write.table(BG, file = "all_facs.tmp.tsv", append = FALSE, quote = FALSE, sep = "\t",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))

idx = grep("wk2",names(BG))
names(BG)[idx]
BG = BG[, -idx]

d1 = BG
names(d1) <- gsub("\\_wk1","", names(d1))
names(d1) <- gsub("X786_O","786_O", names(d1))

x = read.table('facs_exp_design_extended.csv',header=TRUE,sep="\t", as.is = TRUE, quote="", skip = 0, check.names=FALSE);

idx = grep ("wk2", x$sample_name)
x=x[-idx,]

x$sample_name <- gsub("\\_wk1","", x$sample_name)

x$comb_sn =  paste(x$cancer_type2, x$sample_name, sep=".")

names(d1) =  gsub("\\NCI_","", names(d1))

setdiff(names(d1), x$sample_name)

for ( i in 2:(dim(d1)[2])) # add abbreviations CO, BR...etc
{
  idx = which(x$sample_name == names(d1)[i])
  names(d1)[i] = x$comb_sn[idx] 
}

names(d1)

row.names( d1 )  = d1$receptor_esc
d1$receptor_esc = NULL

names(d1) = gsub("\\_","", names(d1))

facs = d1

rownames(facs) <- gsub(" |_|-|\\.|/|,|\\(|)", "", rownames(facs))

rownames_to_remove <- rownames(facs)[grepl("AH|kappa",rownames(facs))]
rownames_to_remove_number <- which(row.names(facs) %in% rownames_to_remove)
facs <- facs[-rownames_to_remove_number,]

## correction of names in facs
names(facs)[names(facs) == "CO.MFIHCT15"] <- "CO.HCT15"
names(facs)[names(facs) == "LC.MFIH23"] <- "LC.H23"
names(facs)[names(facs) == "LE.HL60(TB)"] <- "LE.HL60"
names(facs)[names(facs) == "OV.NCI/ADRRES"] <- "OV.ADRRES"
names(facs)[names(facs) == "RE.786O"] <- "RE.7860"

names(facs)
write.table(facs, "FACS_no_ICS_April2022.tsv", sep = "\t")

####

################################################################

#### transcriptome ####

#td=read.table('transcriptome.tsv',header=TRUE,sep="\t", as.is = TRUE, quote="");
#tmp = td

#a = td[1:2,9:67]
#a=t(a)

 xt = read.table('transcr_exp_design.tsv',header=TRUE,sep="\t", as.is = TRUE, quote="", skip = 0);
 #xt$cancer =  sapply(strsplit(xt$cell_linet,"\\_"), "[[",1)

 xt$clt =  sapply(strsplit(xt$cell_linet,"\\_"), "[[",2)
 
 td2 = read.table('../Curate_ArrayData_April2022/GSE32474_series_matrix_curated.txt',header=TRUE,sep="\t", as.is = TRUE, skip = 0);
td=td2
 
 nm = names(td)
 
 kix =  sapply(strsplit(nm,"\\_"), "[[",1)
 length(unique(kix))
 kix = (intersect(kix, xt$code))
 
 #row.names(td) = td$X
 #td$X = NULL
 
 names(td) =   sapply(strsplit(names(td),"\\_"), "[[",1)
 td =subset(td, select = kix)
 head(td)
 
 for ( i in 1:(dim(td)[2])) # add abbreviations CO, BR...etc
 {
   idx = which(xt$code == names(td)[i])
   names(td)[i] = paste( xt$cancer2[idx], xt$clt[idx], sep=".") 
 }
 
 length(intersect(names(facs), names(td)))
 
 setdiff(names(facs), names(td))

## correction of names in td
 names(td)[names(td) == "BR.MDAMB435"] <- "ME.MDAMB435"
 names(td)[names(td) == "BR.MCF7ADR"] <- "OV.ADRRES"
 
 setdiff(names(facs), names(td)) # only BR.MDAMB468 allowed!

 
 #### only 59 Cell Lines in Transcriptomics data! ####
 
#### ---------------- proteome ####

setwd(wd)

lutp = read.table('prot_nci60_mapping.csv',header=TRUE,sep=",", as.is = TRUE, quote="");
idx = duplicated(lutp$Pnumber)
lutp = lutp[!idx,]

#### NEW FILTERED DATA HERE ####
d1=read.table("../Curate_Px_April2022/proteinGroups_esc_polished_filtered_log10.txt",header=TRUE,sep="\t", as.is = TRUE, check.names = F);
d2 <- d1
intersect(names(facs), names(d2))

length(intersect(names(facs), names(d2)))

setdiff(names(facs), names(d2)) # only "BR.MDAMB468" should appear here!!!

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex(bibstyle="unsrt")


###################################################
### code chunk number 2: load_lib_data
###################################################
library(omicade4)
data(NCI60_4arrays)

names(NCI60_4arrays)
NCI60_4arrays$hgu133p2 <- NCI60_4arrays$hgu133<- NCI60_4arrays$hgu95 <- NULL

NCI60_4arrays$FACS = facs

#NCI60_4arrays$agilent = NULL
  
cmn = intersect(names(facs), names(d2))

idx = grep("LOX", cmn)
cmn = cmn[-idx]

d2 = subset(d2, select = cmn)
# d2[d2==0] = NA
# idx = complete.cases(d2)

#cnts1 <- apply(d2[,1:33], 1, function(x) length(which(x>0)))
#idx = which(cnts1 > 20 )
#d4=d2[idx,]
#d4[is.na(d4)] = 0

#NCI60_4arrays$agilent = subset(NCI60_4arrays$agilent, select = cmn)

NCI60_4arrays$Px = d2

NCI60_4arrays$FACS = subset(NCI60_4arrays$FACS, select = cmn)

NCI60_4arrays$Tx = td

NCI60_4arrays$Tx = subset(NCI60_4arrays$Tx, select = cmn)

NCI60_4arrays$Px = NULL
NCI60_4arrays$FACS = NULL
NCI60_4arrays$Tx = NULL

NCI60_4arrays$Tx = td

NCI60_4arrays$Tx = subset(NCI60_4arrays$Tx, select = cmn)

NCI60_4arrays$Px = d2

NCI60_4arrays$FACS = facs

NCI60_4arrays$FACS = subset(NCI60_4arrays$FACS, select = cmn)

#head(NCI60_4arrays$agilent)


###################################################
### code chunk number 3: check_col_number
###################################################

NCI60_4arrays$agilent <-NULL
sapply(NCI60_4arrays, dim)

facs <- NCI60_4arrays$FACS
library(e1071)
# var(unlist(lapply(facs[grep("PR", colnames(facs))], as.numeric), use.names = FALSE))
cancer_type_uniq <- colnames(NCI60_4arrays$Tx)
cancer_type_uniq <- sapply(strsplit(cancer_type, split="\\."), function(x) x[1])
cancer_type_uniq <- unique(cancer_type_uniq)

for (ct in cancer_type_uniq) {
  print(ct)
  print(kurtosis(unlist(lapply(facs[grep(ct, colnames(facs))], as.numeric), use.names = FALSE)))
  print(skewness(unlist(lapply(facs[grep(ct, colnames(facs))], as.numeric), use.names = FALSE)))
  print(var(unlist(lapply(facs[grep(ct, colnames(facs))], as.numeric), use.names = FALSE)))
  # leptokurtic
  # skewed towards the right
}

kurtosis(unlist(lapply(facs[grep("PR", colnames(facs))], as.numeric), use.names = FALSE))
skewness(unlist(lapply(facs[grep("PR", colnames(facs))], as.numeric), use.names = FALSE))

PR <- facs[grep("PR", colnames(facs))]
PR_ <- rbind(matrix(ncol = 1, nrow = 0), as.matrix(PR[,1]))
PR_ <- rbind(PR_, as.matrix(PR[,2]))
PR_ <- as.data.frame(PR_)
colnames(PR_) <- c("MFI")
PR_$group <- rep("PR", dim(PR_)[1])

BR <- facs[grep("BR", colnames(facs))]
BR_ <- rbind(matrix(ncol = 1, nrow = 0), as.matrix(BR[,1]))
BR_ <- rbind(BR_, as.matrix(BR[,2]))
BR_ <- rbind(BR_, as.matrix(BR[,3]))
BR_ <- rbind(BR_, as.matrix(BR[,4]))
BR_ <- rbind(BR_, as.matrix(BR[,5]))

BR_ <- as.data.frame(BR_)
colnames(BR_) <- c("MFI")
BR_$group <- rep("BR", dim(BR_)[1])

BR_PR_  <- rbind(BR_, PR_)

library(car)
leveneTest(MFI ~ group, data=BR_PR_)

  ###################################################
### code chunk number 4: Check_col_order
###################################################
all(apply((x <- sapply(NCI60_4arrays, colnames))[,-1], 2, function(y)
identical(y, x[,1])))


###################################################
### code chunk number 5: cluster_data_separately
###################################################
library(amap)
library(dendextend)
library(RColorBrewer)
png("dend_spearman.png", height = 1500, width = 1500, pointsize = 15, res = 144)
layout(matrix(1:3, 1, 3))
par(mar=c(7, 1, 0.1, 7))
old_names <- names(NCI60_4arrays)
names(NCI60_4arrays) <- c("Tx", "Px", "FACS")
for (i in 1:length(NCI60_4arrays)) {
  df <- NCI60_4arrays[[i]]
  df_n <- names(NCI60_4arrays[i])
  d <- Dist(t(df), method = "spearman")
  print(d)
  hcl <- hclust(d, method = "ward.D2")
  dend <- as.dendrogram(hcl)
  # https://stackoverflow.com/questions/33683862/first-entry-from-string-split
  colors_to_use <- sapply(strsplit(colnames(df), "\\."), `[`, 1)
  unique_colors <- unique(colors_to_use)
  color_numbers <- brewer.pal(n = 9, name = 'Set1')
  names(color_numbers)=unique_colors
  final_color_vec <- color_numbers[colors_to_use]
  final_color_vec <- final_color_vec[order.dendrogram(dend)]
  labels_colors(dend) <- final_color_vec
  plot(dend, horiz=TRUE, xlab = df_n)
}
dev.off()

svg("dend_spearman.svg", height = 10, width = 10)
layout(matrix(1:3, 1, 3))
par(mar=c(7, 1, 0.1, 7))
for (i in 1:length(NCI60_4arrays)) {
  df <- NCI60_4arrays[[i]]
  df_n <- names(NCI60_4arrays[i])
  d <- Dist(t(df), method = "spearman")
  print(d)
  hcl <- hclust(d, method = "ward.D2")
  dend <- as.dendrogram(hcl)
  # https://stackoverflow.com/questions/33683862/first-entry-from-string-split
  colors_to_use <- sapply(strsplit(colnames(df), "\\."), `[`, 1)
  unique_colors <- unique(colors_to_use)
  color_numbers <- brewer.pal(n = 9, name = 'Set1')
  names(color_numbers)=unique_colors
  final_color_vec <- color_numbers[colors_to_use]
  final_color_vec <- final_color_vec[order.dendrogram(dend)]
  labels_colors(dend) <- final_color_vec
  plot(dend, horiz=TRUE, xlab = df_n)
}
dev.off()

pdf("dend_spearman.pdf", height = 10, width = 10)
layout(matrix(1:3, 1, 3))
par(mar=c(7, 1, 0.1, 7))
for (i in 1:length(NCI60_4arrays)) {
  df <- NCI60_4arrays[[i]]
  df_n <- names(NCI60_4arrays[i])
  d <- Dist(t(df), method = "spearman")
  print(d)
  hcl <- hclust(d, method = "ward.D2")
  dend <- as.dendrogram(hcl)
  # https://stackoverflow.com/questions/33683862/first-entry-from-string-split
  colors_to_use <- sapply(strsplit(colnames(df), "\\."), `[`, 1)
  unique_colors <- unique(colors_to_use)
  color_numbers <- brewer.pal(n = 9, name = 'Set1')
  names(color_numbers)=unique_colors
  final_color_vec <- color_numbers[colors_to_use]
  final_color_vec <- final_color_vec[order.dendrogram(dend)]
  labels_colors(dend) <- final_color_vec
  plot(dend, horiz=TRUE, xlab = df_n)
}
dev.off()

# png("dend_euclidean.png", height = 1500, width = 2000, pointsize = 15, res = 144)
# layout(matrix(1:3, 1, 3))
# par(mar=c(7, 1, 0.1, 7))
# for (i in 1:length(NCI60_4arrays)) {
#   df <- NCI60_4arrays[[i]]
#   df_n <- names(NCI60_4arrays[i])
#   d <- dist(t(df), method = "euclidean")
#   print(d)
#   hcl <- hclust(d, method = "ward.D2")
#   dend <- as.dendrogram(hcl)
#   # https://stackoverflow.com/questions/33683862/first-entry-from-string-split
#   colors_to_use <- sapply(strsplit(colnames(df), "\\."), `[`, 1)
#   unique_colors <- unique(colors_to_use)
#   color_numbers <- brewer.pal(n = 9, name = 'Set1')
#   names(color_numbers)=unique_colors
#   final_color_vec <- color_numbers[colors_to_use]
#   final_color_vec <- final_color_vec[order.dendrogram(dend)]
#   labels_colors(dend) <- final_color_vec
#   plot(dend, horiz=TRUE, xlab = df_n)
# }
# dev.off()

#library(bioDist)
 
#png("dend_spearman_bioDist.png", height = 1500, width = 1500, pointsize = 15, res = 144)
#layout(matrix(1:3, 1, 3))
#par(mar=c(7, 1, 0.1, 7))
#for (i in 1:length(NCI60_4arrays)) {
#  df <- NCI60_4arrays[[i]]
#  df_n <- names(NCI60_4arrays[i])
#  d <- spearman.dist(t(df))
#  print(d)
#  hcl <- hclust(d, method = "ward.D2")
#  dend <- as.dendrogram(hcl)
#  plot(dend, horiz=TRUE, xlab = df_n)
#}
#dev.off()

# png("dend_eclidean_avg_link.png", height = 1500, width = 2000, pointsize = 15, res = 144)
# layout(matrix(1:3, 1, 3))
# par(mar=c(7, 1, 0.1, 7))
# for (i in 1:length(NCI60_4arrays)) {
#   df <- NCI60_4arrays[[i]]
#   df_n <- names(NCI60_4arrays[i])
#   d <- dist(t(df), method = "euclidean")
#   print(d)
#   hcl <- hclust(d, method = "average")
#   dend <- as.dendrogram(hcl)
#   # https://stackoverflow.com/questions/33683862/first-entry-from-string-split
#   colors_to_use <- sapply(strsplit(colnames(df), "\\."), `[`, 1)
#   unique_colors <- unique(colors_to_use)
#   color_numbers <- brewer.pal(n = 9, name = 'Set1')
#   names(color_numbers)=unique_colors
#   final_color_vec <- color_numbers[colors_to_use]
#   final_color_vec <- final_color_vec[order.dendrogram(dend)]
#   labels_colors(dend) <- final_color_vec
#   plot(dend, horiz=TRUE, xlab = df_n)
# }
# dev.off()
#names(NCI60_4arrays) <- old_names

###################################################
### code chunk number 6: MCIA
###################################################
mcoin <- mcia(NCI60_4arrays, cia.nf=10)
class(mcoin)

###################################################
### code chunk number 7: get_cancertype
###################################################
cancer_type <- colnames(NCI60_4arrays$Tx)
cancer_type <- sapply(strsplit(cancer_type, split="\\."), function(x) x[1])
cancer_type


###################################################
### code chunk number 8: plot_mcia_cancertype
###################################################
library(lattice)
source("my_omicade4.R")
# trellis.device(device="png", filename="no_lox_1.png", width=1500, height=1500)
png("mcia.png", height = 1500, width = 1500, pointsize = 15, res = 144)
# colcode <- sapply(strsplit(colnames(NCI60_4arrays$array), split="\\."), 
#                   function(x) x[1])
# x11()
# plot(mcoin, sample.lab=T, sample.color=as.factor(colcode))

my.plot.mcia(mcoin, axes=1:2, phenovec=cancer_type, sample.lab=F, df.color=1:3)

dev.off()

svg("mcia.svg",width=10,height=10)
my.plot.mcia(mcoin, axes=1:2, phenovec=cancer_type, sample.lab=FALSE, df.color=1:3)

dev.off()

pdf("mcia.pdf",width=10,height=10)
my.plot.mcia(mcoin, axes=1:2, phenovec=cancer_type, sample.lab=FALSE, df.color=1:3)

dev.off()

## calculate 
# and blue line corresponds to the percentage of variance of each PC, calculated as the eigenvalue divided by sum of all eigenvalues.
sum(mcoin$mcoa$pseudoeig[1:10])/sum(mcoin$mcoa$pseudoeig) # 0.7471196


###################################################
### code chunk number 9: selectvar_mcia_melan
###################################################
#CD105 <- selectVar(mcoin, a1.lim=c(-Inf, 0), a2.lim=c(2, Inf))
# melan_gene <- selectVar(mcoin, a1.lim=c(1, Inf))
#CD105

#colon <- selectVar(mcoin, a1.lim=c(0.5, Inf), a2.lim = c(-Inf, -0.5))

#head(colon)
###################################################
### code chunk number 10: plot_mcia_axes
###################################################
plot(mcoin, axes=c(1, 2), phenovec=cancer_type, sample.lab=F, df.color=1:3, sample.legend = T)

png("mcia_axes_1,2_June2021.png", height = 1500, width = 2000, pointsize = 15, res = 144)
my.plot.mcia(mcoin, axes=c(1, 2), phenovec=cancer_type, sample.lab=T, df.color=1:3)
dev.off()
plot(mcoin, axes=c(1, 3), phenovec=cancer_type, sample.lab=T, df.color=1:3)

png("mcia_axes_1,3_June2021.png", height = 1500, width = 2000, pointsize = 15, res = 144)
my.plot.mcia(mcoin, axes=c(1, 3), phenovec=cancer_type, sample.lab=T, df.color=1:3)
dev.off()

plot(mcoin, axes=c(2, 3), phenovec=cancer_type, sample.lab=FALSE, df.color=1:3)
plot(mcoin, axes=c(1, 4), phenovec=cancer_type, sample.lab=FALSE, df.color=1:3, )


###################################################
### code chunk number 11: plotvar_mcia
###################################################
# trellis.device(device="png", filename="no_lox_scores.png", width=540, height=540)
#png("no_lox_scores.png", height = 1500, width = 1500, pointsize = 15, res = 144)
#geneStat <- plotVar(mcoin, var=c("CD100", "CCR10"), var.lab=TRUE)
#dev.off()

## S3 method for class 'mcia'
top5_axis1 <- topVar(mcoin, axis = 1, end = "both", topN = 5)
write.table(top5_axis1, file = "facs_schindler_multiomics_top10_hits_x_axis.tsv", append = FALSE, quote = FALSE, sep = "\t",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))
top5_axis2 <- topVar(mcoin, axis = 2, end = "both", topN = 5)
write.table(top5_axis2, file = "facs_schindler_multiomics_top10_hits_y_axis.tsv", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))
top10 <- c(as.character(unlist(top5_axis1)), as.character(unlist(top5_axis2)))
png("facs_schindler_multiomics_mcia_top10_filtered.png", height = 1500, width = 1500, pointsize = 15, res = 144)
plotVar(mcoin, var=top10, var.lab=TRUE, sepID.data=4, sepID.sep="\\.")
dev.off()
svg("facs_schindler_multiomics_mcia_top10_filtered.svg", height = 10, width = 10)
plotVar(mcoin, var=top10, var.lab=TRUE, sepID.data=4, sepID.sep="\\.")
dev.off()

top10_axis1 <- topVar(mcoin, axis = 1, end = "both", topN = 10)
write.table(top10_axis1, file = "facs_schindler_multiomics_top20_hits_x_axis.tsv", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))
top10_axis2 <- topVar(mcoin, axis = 2, end = "both", topN = 10)
write.table(top10_axis2, file = "facs_schindler_multiomics_top20_hits_y_axis.tsv", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))
top20 <- c(as.character(unlist(top10_axis1)), as.character(unlist(top10_axis2)))
png("facs_schindler_multiomics_mcia_top20_filtered.png", height = 1500, width = 1500, pointsize = 15, res = 144)
plotVar(mcoin, var=top20, var.lab=TRUE, sepID.data=4, sepID.sep="\\.")
dev.off()
svg("facs_schindler_multiomics_mcia_top20_filtered.svg", height = 10, width = 10)
plotVar(mcoin, var=top20, var.lab=TRUE, sepID.data=4, sepID.sep="\\.")
dev.off()

top20_axis1_positive <- topVar(mcoin, axis = 1, end = "positive", topN = 20)
write.table(top10_axis1, file = "facs_schindler_multiomics_top20_hits_x_axis_positive.tsv", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))
png("facs_schindler_multiomics_mcia_top20_filtered_positive.png", height = 1500, width = 1500, pointsize = 15, res = 144)
plotVar(mcoin, var=as.character(unlist(top20_axis1_positive)), var.lab=TRUE, sepID.data=4, sepID.sep="\\.")
dev.off()
svg("facs_schindler_multiomics_mcia_top20_filtered_positive.svg", height = 10, width = 10)
plotVar(mcoin, var=as.character(unlist(top20_axis1_positive)), var.lab=TRUE, sepID.data=4, sepID.sep="\\.")
dev.off()

## CUT OUT OUR HITS

colon <- selectVar(mcoin, a1.lim=c(-1, -0.4), a2.lim=c(-1, -0.38))
colon <- colon[colon$FACS,]
write.table(colon, file = "facs_schindler_multiomics_colon_hits.tsv", append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))
png("facs_schindler_multiomics_mcia_colon.png", height = 1500, width = 1500, pointsize = 15, res = 144)
plotVar(mcoin, var=as.character(unlist(colon)), var.lab=F)
dev.off()
svg("facs_schindler_multiomics_mcia_colon.svg", height = 10, width = 10)
plotVar(mcoin, var=as.character(unlist(colon)), var.lab=F)
dev.off()

leuk <- selectVar(mcoin, a1.lim=c(-Inf, -1), a2.lim=c(0, Inf))
leuk <- leuk[leuk$FACS,]
write.table(leuk, file = "facs_schindler_multiomics_leuk_hits.tsv", append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))
png("facs_schindler_multiomics_mcia_leuk.png", height = 1500, width = 1500, pointsize = 15, res = 144)
plotVar(mcoin, var=as.character(unlist(leuk)), var.lab=F)
dev.off()
svg("facs_schindler_multiomics_mcia_leuk.svg", height = 10, width = 10)
plotVar(mcoin, var=as.character(unlist(leuk)), var.lab=F)
dev.off()

me <- selectVar(mcoin, a1.lim=c(-0.25, 0.15), a2.lim=c(0.5, Inf))
me <- me[me$FACS,]
write.table(me, file = "facs_schindler_multiomics_me_hits.tsv", append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))
png("facs_schindler_multiomics_mcia_me.png", height = 1500, width = 1500, pointsize = 15, res = 144)
plotVar(mcoin, var=as.character(unlist(me)), var.lab=F)
dev.off()
svg("facs_schindler_multiomics_mcia_me.svg", height = 10, width = 10)
plotVar(mcoin, var=as.character(unlist(me)), var.lab=F)
dev.off()

cns <- selectVar(mcoin, a1.lim=c(0.3, Inf), a2.lim=c(0.5, Inf))
cns <- cns[cns$FACS,]
write.table(cns, file = "facs_schindler_multiomics_cns_hits.tsv", append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))
png("facs_schindler_multiomics_mcia_cns.png", height = 1500, width = 1500, pointsize = 15, res = 144)
plotVar(mcoin, var=as.character(unlist(cns)), var.lab=F)
dev.off()
svg("facs_schindler_multiomics_mcia_cns.svg", height = 10, width = 10)
plotVar(mcoin, var=as.character(unlist(cns)), var.lab=F)
dev.off()

re <- selectVar(mcoin, a1.lim=c(0.25, Inf), a2.lim=c(-Inf, -0.5))
re <- re[re$FACS,]
write.table(re, file = "facs_schindler_multiomics_re_hits.tsv", append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))
png("facs_schindler_multiomics_mcia_re.png", height = 1500, width = 1500, pointsize = 15, res = 144)
plotVar(mcoin, var=as.character(unlist(re)), var.lab=F)
dev.off()
svg("facs_schindler_multiomics_mcia_re.svg", height = 10, width = 10)
plotVar(mcoin, var=as.character(unlist(re)), var.lab=F)
dev.off()

harmonized_ids <- read.table("../Homogenize_Identifiers/harmonized_final.ids.csv", header = T, sep = "\t")
harmonized_ids$receptor_esc <- gsub(" |_|-|\\.|/|,|\\(|)", "", harmonized_ids$receptor_esc)

leuk_gene <- merge(harmonized_ids, leuk, by.x = "receptor_esc", by.y = "var")
leuk_gene <- leuk_gene[, names(leuk_gene) %in% c("receptor_esc", "approved_receptor")]
names(leuk_gene) <- c("receptor", "gene_symbol")
write.table(leuk_gene, file = "leuk_hits_gene_symbol.tsv", append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))

me_gene <- merge(harmonized_ids, me, by.x = "receptor_esc", by.y = "var")
me_gene <- me_gene[, names(me_gene) %in% c("receptor_esc", "approved_receptor")]
names(me_gene) <- c("receptor", "gene_symbol")
write.table(me_gene, file = "me_hits_gene_symbol.tsv", append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))

colon_gene <- merge(harmonized_ids, colon, by.x = "receptor_esc", by.y = "var")
colon_gene <- colon_gene[, names(colon_gene) %in% c("receptor_esc", "approved_receptor")]
names(colon_gene) <- c("receptor", "gene_symbol")
write.table(colon_gene, file = "colon_hits_gene_symbol.tsv", append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))

re_gene <- merge(harmonized_ids, re, by.x = "receptor_esc", by.y = "var")
re_gene <- re_gene[, names(re_gene) %in% c("receptor_esc", "approved_receptor")]
names(re_gene) <- c("receptor", "gene_symbol")
write.table(re_gene, file = "re_hits_gene_symbol.tsv", append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))

cns_gene <- merge(harmonized_ids, cns, by.x = "receptor_esc", by.y = "var")
cns_gene <- cns_gene[, names(cns_gene) %in% c("receptor_esc", "approved_receptor")]
names(cns_gene) <- c("receptor", "gene_symbol")
write.table(cns_gene, file = "cns_hits_gene_symbol.tsv", append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"))
save.image("mcia_April2022_final.Rdata")
