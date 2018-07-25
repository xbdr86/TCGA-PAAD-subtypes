setwd("~/Desktop/Gut-Jeroni/")
require(gdata)
dyrk1a <- read.xls ("dyrk1a.xlsx", sheet = "dyrk1a", header = TRUE)
metadata <- read.xls ("dyrk1a.xlsx", sheet = "metadata", header = TRUE)
dyrk1a <- merge(dyrk1a, metadata, by="Title")
tumor <- dyrk1a[which(dyrk1a$Group=="tumor"),]
adjacent <- dyrk1a[which(dyrk1a$Group=="adj"),]
write.table(tumor, "tumor.txt", sep="\t", append = FALSE)
write.table(adjacent, "adjacent.txt", sep="\t", append = FALSE)
