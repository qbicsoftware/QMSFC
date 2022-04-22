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

idx = grep("wk2|wk1",names(BG))
names(BG)[idx]
BG = BG[, idx]

names(BG)[names(BG) == "786_O_wk1"] <- "X786_wk1"
names(BG)[names(BG) == "786_O_wk2"] <- "X786_wk2"

rownames(BG) <- d1$receptor_esc

library("ggpubr")

output <- matrix(ncol = 4, nrow = 0)
even_cols <- seq(2, ncol(BG), by = 2)
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r

for (even_col in even_cols) {
  # test both cols with shapiro and scream if they should be normally distributed!
  shap1 <- shapiro.test(BG[,even_col])
  stopifnot(shap1$p.value < 0.01)
  shap2 <- shapiro.test(BG[,even_col - 1])
  stopifnot(shap2$p.value < 0.01)
  cell_line <- colnames(BG)[even_col - 1]
  #cell_line <- gsub("\\_wk1","", cell_line)
  cell_line <- gsub("MFI\\_","", cell_line)
  if (cell_line == "786_O_wk1") {
    print("SHFLSDJHK")
    cell_line <- "786_wk1"
  }
  pdf(paste(cell_line, "corr_spearman_wk1_vs_wk2.pdf", sep = "_"))
  print(ggscatter(BG, x = colnames(BG)[even_col - 1], y = colnames(BG)[even_col], 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "wk1", ylab = "wk2", main = paste("Correllation test of ", cell_line, " week 1 versus week 2")))
  dev.off()
  # https://statisticsbyjim.com/basics/correlations/
  res <-cor.test(BG[,even_col - 1], BG[,even_col],  method = "spearman", exact = F)
  log2_FC <- log2(BG[,even_col]) - log2(BG[,even_col - 1])
  d2 <- as.data.frame(log2_FC)
  rownames(d2) <- rownames(BG)
  colnames(d2) <- c(paste(cell_line, "_wk1_vs_wk2_log2FC", sep = ""))
  write.table(d2, paste(cell_line, "_wk1_vs_wk2_log2FC.tsv", sep = ""), sep = "\t")
  # table with cell line name, r-value, p-value, number of log2FCs greater |2|
  output <- rbind(output, c(cell_line, as.numeric(res$estimate), res$p.value, sum(log2_FC < -2 | log2_FC > 2)))
  #print(sum(log2_FC < -0.5 | log2_FC > 0.5))
  #print(sum(log2_FC < -2 | log2_FC > 2))
  #print("#")
  #Sys.sleep(2)
}

output <- as.data.frame(output)
colnames(output) <- c("cell_line", "rho", "p-value", "number_of_|log2FC|>2")
write.table(output, "wk1_vs_wk2_correlation_log2FC_analysis.tsv",sep = "\t", row.names = F)
