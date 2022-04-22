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

facs_ids <- read.table("NCI60_leukemia_Clean.tsv", header=TRUE ,sep="\t", as.is = TRUE, quote="", skip = 0, check.names=FALSE)
recs_clean <- c()
for (i in 2:length(facs_ids$receptor_esc)) {
  rec_esc <- facs_ids$receptor_esc[i]
  rec_esc_split <- strsplit(rec_esc, " ")
  rec_clean <- rec_esc_split[[1]][1]
  recs_clean <- c(recs_clean, rec_clean)
}

cd_ids <- read.table("gene_group-CD_molecules.tsv", header=TRUE ,sep="\t", as.is = TRUE, quote="", skip = 1, check.names=FALSE)
facs_clean <- facs_ids
facs_clean$receptor <- c(facs_ids$receptor[1], recs_clean)
facs_clean <- facs_clean[,colnames(facs_clean) == "receptor" | colnames(facs_clean) == "receptor_esc"]
colnames(facs_clean) <- c("receptor", "approved_receptor")
facs_clean$approved_receptor <- rep("NF", length(facs_clean))

for (i in 2:length(facs_clean$receptor)) {
  appr_grep <- which(facs_clean$receptor[i] == cd_ids$`Approved symbol`)
  if (length(appr_grep) > 0) {
    facs_clean[i, 2] <- cd_ids$`Approved symbol`[appr_grep]
  } else {
    for (k in 1:length(cd_ids$Synonyms)) {
      synonyms_split <- strsplit(cd_ids$Synonyms[k], ",")[[1]]
      if (length(synonyms_split) > 0) {
        for (j in 1:length(synonyms_split)) {
          is_match <- which(facs_clean$receptor[i] == synonyms_split[j])
          if (length(is_match) > 0) {
            facs_clean[i, 2] <- cd_ids$`Approved symbol`[k]
          }
        }
      }
    }
  }
}

facs_clean$receptor_esc <- facs_ids$receptor_esc
write.table(facs_clean, file = "pre.harmonized.ids.tsv", sep = "\t", col.names = T, row.names = F)
