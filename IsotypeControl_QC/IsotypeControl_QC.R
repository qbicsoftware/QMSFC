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

isotypes <- read.table("../Homogenize_Identifiers/isotype.ids.csv", as.is = T, sep = "\t", stringsAsFactors = F, header = T)

#### --------------FACS ####
wd <- getwd()
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
facs_iso <- facs[rownames_to_remove_number,]
facs <- facs[-rownames_to_remove_number,]
facs_clean <- facs

#idx = grep("LOX", names(facs_clean))
#facs_clean = facs_clean[,-idx]

#ms_igg3 <- as.numeric(facs_iso[facs_iso$srt == "297", -c(1:4)])

#cd167 <- facs_clean[rownames(facs_clean) == "CD167a (DDR1) ",]

#cd167_iso <- cd167 - ms_igg3

#SF <- facs_clean[,c(1:2)]

#### boxplot and histogram by isotype ####

library(RColorBrewer)
colors_to_use <- sapply(strsplit(colnames(facs_clean), "\\."), `[`, 1)
unique_colors <- unique(colors_to_use)
color_numbers <- brewer.pal(n = 9, name = 'Set1')
names(color_numbers)=unique_colors
final_color_vec <- color_numbers[colors_to_use]

# for (i in 1:length(row.names(facs_iso))) {
#   rec <- row.names(facs_iso)[i]
#   print(rec)
#   png(paste("Boxplot_MFI_IsotypeControl_", rec, ".png", collapse = ""), height = 1500, width = 2750, units = "px", pointsize = 30)
#   par(mar = c(10,4,2,2))
#   boxplot(log10(facs_iso[i,]), main = paste("MFI values of IsotypeControl ", rec, ".png", collapse = ""), las = 2, ylab = "MFI")
#   dev.off()
#   png(paste("Histogram_MFI_IsotypeControl_", rec, ".png", collapse = ""), height = 1500, width = 2750, units = "px", pointsize = 30)
#   hist(as.numeric(facs_iso[i,]), main = paste("MFI values of IsotypeControl ", rec, collapse = ""), 
#        breaks = 200, xlab = "MFI")
#   dev.off()
# }

#### boxplot and histogram by isotype mouse ####
facs_iso_mm <- facs_iso[grep("Ms",row.names(facs_iso)),]
png("Boxplot_MFI_IsotypeControl_Mouse.png", height = 1500, width = 2750, units = "px", pointsize = 30)
par(mar = c(10,4,2,2))
boxplot(facs_iso_mm, main = "MFI values of IsotypeControl Mouse Only", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()
pdf("Boxplot_MFI_IsotypeControl_Mouse.pdf", width = 12)
par(mar = c(10,4,2,2))
boxplot(facs_iso_mm, main = "MFI values of IsotypeControl Mouse Only", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()
#png("Histogram_MFI_IsotypeControl_Mouse.png", height = 1500, width = 2750, units = "px", pointsize = 30)
#hist(as.numeric(unlist(facs_iso_mm)), main = "MFI values of IsotypeControl Mouse Only", 
#     breaks = 200, xlab = "MFI")
#dev.off()

#### boxplot of MsIgG3kappaITCL ####
facs_iso_msigg3 <- facs_iso[grep("MsIgG3",row.names(facs_iso)),]
png("Boxplot_MFI_IsotypeControl_Mouse_MsIgG3.png", height = 1500, width = 2750, units = "px", pointsize = 30)
par(mar = c(10,4,2,2))
boxplot(facs_iso_msigg3, main = "MFI values of IsotypeControl Mouse MsIgG3kappaITCL Only", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()
pdf("Boxplot_MFI_IsotypeControl_Mouse_MsIgG3.pdf", width = 12)
par(mar = c(10,4,2,2))
boxplot(facs_iso_msigg3, main = "MFI values of IsotypeControl Mouse MsIgG3kappaITCL Only", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()

#### boxplot of Mouse without MsIgG3kappaITCL
facs_iso_mouse_no_msigg3 <- facs_iso[-grep("MsIgG3",row.names(facs_iso)),]
facs_iso_mouse_no_msigg3 <- facs_iso_mouse_no_msigg3[grep("Ms", row.names(facs_iso_mouse_no_msigg3)),]
png("Boxplot_MFI_IsotypeControl_Mouse_No_MsIgG3.png", height = 1500, width = 2750, units = "px", pointsize = 30)
par(mar = c(10,4,2,2))
boxplot(facs_iso_mouse_no_msigg3, main = "MFI values of IsotypeControl Mouse without MsIgG3kappaITCL", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()
pdf("Boxplot_MFI_IsotypeControl_Mouse_No_MsIgG3.pdf", width = 12)
par(mar = c(10,4,2,2))
boxplot(facs_iso_mouse_no_msigg3, main = "MFI values of IsotypeControl Mouse without MsIgG3kappaITCL", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()

#### boxplot and histogram by isotype rat ####
facs_iso_rat <- facs_iso[grep("Rat",row.names(facs_iso)),]
png("Boxplot_MFI_IsotypeControl_Rat.png", height = 1500, width = 2750, units = "px", pointsize = 30)
par(mar = c(10,4,2,2))
boxplot(facs_iso_rat, main = "MFI values of IsotypeControl Rat Only", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()
pdf("Boxplot_MFI_IsotypeControl_Rat.pdf", width = 12)
par(mar = c(10,4,2,2))
boxplot(facs_iso_rat, main = "MFI values of IsotypeControl Rat Only", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()
#png("Histogram_MFI_IsotypeControl_Rat.png", height = 1500, width = 2750, units = "px", pointsize = 30)
#hist(as.numeric(unlist(facs_iso_rat)), main = "MFI values of IsotypeControl Rat Only", 
#     breaks = 200, xlab = "MFI")
#dev.off()

#### boxplot of isotype AHIg ####
facs_iso_ah <- facs_iso[grep("AH",row.names(facs_iso)),]
png("Boxplot_MFI_IsotypeControl_AHIgGITCL.png", height = 1500, width = 2750, units = "px", pointsize = 30)
par(mar = c(10,4,2,2))
boxplot(facs_iso_ah, main = "MFI values of IsotypeControl AHIgGITCL Only", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()
pdf("Boxplot_MFI_IsotypeControl_AHIgGITCL.pdf", width = 12)
par(mar = c(10,4,2,2))
boxplot(facs_iso_ah, main = "MFI values of IsotypeControl AHIgGITCL Only", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()

#### boxplot and histogram of all isotypes #### 
png("Boxplot_MFI_IsotypeControl_All.png", height = 1500, width = 2750, units = "px", pointsize = 30)
par(mar = c(10,4,2,2))
boxplot(facs_iso, main = "MFI values of IsotypeControl All", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()
pdf("Boxplot_MFI_IsotypeControl_All.pdf", width = 12)
par(mar = c(10,4,2,2))
boxplot(facs_iso, main = "MFI values of IsotypeControl All", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()
# png("Histogram_MFI_IsotypeControl_All.png", height = 1500, width = 2750, units = "px", pointsize = 30)
# hist(as.numeric(unlist(facs_iso[,-c(1,2,3,4)])), main = "MFI values of IsotypeControl All", 
#      breaks = 200, xlab = "MFI")
# dev.off()

#### boxplot of all isotypes except Mouse igG3
facs_iso_no_msigg3 <- facs_iso[-grep("MsIgG3",row.names(facs_iso)),]
png("Boxplot_MFI_IsotypeControl_All_No_MsIgG3.png", height = 1500, width = 2750, units = "px", pointsize = 30)
par(mar = c(10,4,2,2))
boxplot(facs_iso_no_msigg3, main = "MFI values of IsotypeControl All without MsIgG3kappaITCL", las = 2, ylab = "MFI", col = final_color_vec)
dev.off()
pdf("Boxplot_MFI_IsotypeControl_All_No_MsIgG3.pdf", width = 12)
par(mar = c(10,4,2,2))
boxplot(facs_iso_no_msigg3, main = "MFI values of IsotypeControl All without MsIgG3kappaITCL", las = 2, ylab = "MFI", col = final_color_vec)
#abline(h = 1)
#abline(h = 2)
dev.off()

sum(facs_iso_no_msigg3 <= 1)
#[1] 368
sum(facs_iso_no_msigg3 > 1)
#[1] 172
sum(facs_iso_no_msigg3 > 2)
#[1] 62
sum(facs_iso_no_msigg3 <= 2)
#[1] 478

#### boxplot and histogram of all facs MFI #### 
# png("Boxplot_MFI_FACS_All.png", height = 1500, width = 2750, units = "px", pointsize = 30)
# par(mar = c(10,4,2,2))
# boxplot(facs_clean, main = "MFI values of FACS All", las = 2, ylab = "MFI")
# dev.off()

png("Boxplot_MFI_FACS_All_log10.png", height = 1500, width = 2750, units = "px", pointsize = 30)
par(mar = c(10,4,2,2))
boxplot(log10(facs_clean), main = "MFI values of FACS All no ICs", las = 2, ylab = "MFI log10")
dev.off()
pdf("Boxplot_MFI_FACS_All_log10.pdf", width = 12)
par(mar = c(10,4,2,2))
boxplot(log10(facs_clean), main = "MFI values of FACS All no ICs", las = 2, ylab = "MFI log10", outcex=0.25, col = final_color_vec)
dev.off()
# png("Histogram_MFI_FACS_All.png", height = 1500, width = 2750, units = "px", pointsize = 30)
# hist(as.numeric(unlist(facs_clean)), main = "MFI values of FACS All", 
#      breaks = 200, xlab = "MFI")
# dev.off()

