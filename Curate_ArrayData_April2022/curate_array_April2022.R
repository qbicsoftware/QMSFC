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

require("tidyverse")

array <- read.table('GSE32474_series_matrix.txt',header=TRUE,sep="\t", as.is = TRUE, quote="", skip = 75, check.names = F)

colnames(array) <- str_remove_all(colnames(array), "\"")
array[,1] <- str_remove_all(array[,1], "\"")
rownames(array) <- array[,1]
array <- array[,-1]

id_mapping <- read.table('martquery_0226141000_361.txt', header = T, sep = "\t", as.is = T, quote = "", check.names = F)
# remove probes where no HGNC symbol exists
id_mapping <- id_mapping[!id_mapping$`HGNC symbol` == "",]
unique_probes <- unique(id_mapping$`AFFY HG U133 Plus 2 probe`)

hgnc2affy <- data.frame(HGNC_Symbol=character(),
                 AFFY_HG_U133_Plus2_Probe=character(), 
                 stringsAsFactors=FALSE) 
# 
for (i in 1:length(unique_probes)) {
  unique_probe <- unique_probes[i]
  hgnc <- c()
  print(i)
  idx <- which(id_mapping$`AFFY HG U133 Plus 2 probe` == unique_probe)
  for (j in 1:length(idx)) {
    affy_hgnc <- id_mapping[idx[j],6]
    if (!affy_hgnc %in% hgnc)
    hgnc <- c(hgnc, affy_hgnc)
  }
  if (length(hgnc) == 1) {
    # https://stackoverflow.com/questions/28467068/add-row-to-dataframe
    # hgnc_str <- paste(hgnc, collapse = " /// ")
    # hgnc2affy[nrow(hgnc2affy) + 1,] = list(hgnc_str, unique_probe)
    hgnc2affy[nrow(hgnc2affy) + 1,] = list(hgnc, unique_probe)
  }
}

## number of rows hgnc2affy without distinct ids
#39984
## number of rows hgnc2affy with distinct ids
#37989
## number of unique symbols without distinct ids
#20974
## number of unique symbols with distinct ids
#19473

save.image("./curate_array.RData")

# Now combine via the exact same HGNC and summarize via mean
unique_symbols <- unique(hgnc2affy$HGNC_Symbol)
array_summarized <- array[1,]
array_summarized <- array_summarized[-1,]
collected_rownames <- c()
for (i in 1:length(unique_symbols)) {
  print(i)
  unique_symbol <- unique_symbols[i]
  idx_hgnc2affy <- which(unique_symbol == hgnc2affy$HGNC_Symbol)
  probes <- hgnc2affy$AFFY_HG_U133_Plus2_Probe[idx_hgnc2affy]
  # collect indices of probes
  idx <- c()
  for (j in 1:length(probes)) {
    probe <- probes[j]
    idx <- c(idx, which(probe == rownames(array)))
  }
  array_summarized[i,] <- colMeans(array[idx,])
  collected_rownames <- c(collected_rownames, unique_symbols[i])
}
rownames(array_summarized) <- collected_rownames
save.image("./curate_array.RData")

write.table(x = array_summarized, "GSE32474_series_matrix_curated.txt", col.names = T, row.names = T, sep = "\t")

