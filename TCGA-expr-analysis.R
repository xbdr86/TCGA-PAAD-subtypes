setwd("~/Desktop/normalized gene expression TCGA/")
metadata <- read.table("all-tumor-manifest.csv", header = TRUE, sep = ",")

setwd("~/Desktop/normalized gene expression TCGA/miRNA_expression/")
miRNA_paad <- read.table("PAAD.normalized.miRNA.tsv", header = TRUE, sep = "\t")

setwd("~/Desktop/normalized gene expression TCGA/gene_expression/")
RNA_paad <- read.table("PAAD.normalized.FPKM-UQ.tsv", header = TRUE, sep = "\t")
dyrk1a <- 


setwd("~/Desktop/normalized gene expression TCGA/")
merged_paad <- read.table("PAAD.normalized.tsv", header = TRUE, sep = "\t")


linearMod <- lm(log(merged_paad$ENSG00000157540.18) ~ log(merged_paad$hsa.mir.148a), data=merged_paad)
plot(log(merged_paad$hsa.mir.148a), log(merged_paad$ENSG00000157540.18))
abline(linearMod)

mean(rna$ENSG00000105976.13) #c-met
mean(rna$ENSG00000157540.18) #dyrk1a
mean(rna$ENSG00000105204.12) #dyrk1b

normal <- metadata[which(metadata$DISEASE.ABBV=="PAAD" & metadata$SAMPLE.TYPE=="SOLID.TISSUE.NORMAL"),]
tumoral <- metadata[which(metadata$DISEASE.ABBV=="PAAD" & metadata$SAMPLE.TYPE=="PRIMARY.TUMOR"),]

rna <- merged_paad[complete.cases(merged_paad), ]
plot(log(rna$ENSG00000157540.18), log(rna$ENSG00000105204.12), xlab = "log(DYRK1A expression)", ylab ="log(DYRK1B expression)")
linearMod <- lm(log(rna$ENSG00000157540.18) ~ log(rna$ENSG00000105204.12), data=rna)
abline(linearMod)
plot(log(rna$ENSG00000157540.18), log(rna$ENSG00000105976.13), xlab = "log(DYRK1A expression)", ylab ="log(c-MET expression)")
linearMod <- lm(log(rna$ENSG00000157540.18) ~ log(rna$ENSG00000105976.13), data=rna)
abline(linearMod)

setwd("~/Desktop/Gut-Jeroni/")
require(gdata)
subtypes <- read.xls ("PDAC-subtypes.xlsx", sheet = "FreezeSamples", header = TRUE)
subtypes$sample_id <- subtypes$Tumor.Sample.ID
subtypes <- merge(subtypes, rna, by="sample_id")
subtypes <- subtypes[order(subtypes$ENSG00000157540.18),] 
moffit.basal <- subtypes[which(subtypes$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX=="1"),]
moffit.classic <- subtypes[which(subtypes$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX=="2"),]
bailey.progenitor <- subtypes[which(subtypes$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX=="3"),]
bailey.adex <- subtypes[which(subtypes$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX=="4"),]
mean(moffit.basal$ENSG00000157540.18) #dyrk1a
sd(moffit.basal$ENSG00000157540.18)
mean(moffit.classic$ENSG00000157540.18) #dyrk1a
sd(moffit.classic$ENSG00000157540.18)
mean(bailey.progenitor$ENSG00000157540.18) #dyrk1a
sd(bailey.progenitor$ENSG00000157540.18)
mean(bailey.adex$ENSG00000157540.18) #dyrk1a
sd(bailey.adex$ENSG00000157540.18)

mean(moffit.basal$Days.to.death, na.rm = TRUE) #dyrk1a
sd(moffit.basal$Days.to.death, na.rm = TRUE)
mean(moffit.classic$Days.to.death, na.rm = TRUE) #dyrk1a
sd(moffit.classic$Days.to.death, na.rm = TRUE)
mean(bailey.progenitor$Days.to.death, na.rm = TRUE) #dyrk1a
sd(bailey.progenitor$Days.to.death, na.rm = TRUE)
mean(bailey.adex$Days.to.death, na.rm = TRUE) #dyrk1a
sd(bailey.adex$Days.to.death, na.rm = TRUE)

mean(moffit.classic[c(1:13),32], na.rm = TRUE) #dyrk1a
sd(moffit.classic[c(1:13),32], na.rm = TRUE) #dyrk1a
mean(moffit.classic[c(15:27),32], na.rm = TRUE) #dyrk1a
sd(moffit.classic[c(15:27),32], na.rm = TRUE) #dyrk1a

aggdata <-aggregate(subtypes$ENSG00000157540.18, by=list(subtypes$Number.of.lymph.nodes), FUN=mean, na.rm=TRUE)
aggdata2 <-aggregate(subtypes$ENSG00000157540.18, by=list(subtypes$Number.of.lymph.nodes), FUN=sd, na.rm=TRUE)
aggdata3 <-aggregate(subtypes$ENSG00000157540.18, by=list(subtypes$Number.of.lymph.nodes), FUN = function(x) c(n = length(x)))
aggregate.total <- merge(aggdata, aggdata2, by="Group.1")
aggregate.total <- merge(aggregate.total, aggdata3, by="Group.1")



rna$SAMPLE.ID <- rna$sample_id
normal <-merge(rna, normal, by="SAMPLE.ID")
mean(normal$ENSG00000105976.13) #c-met
mean(normal$ENSG00000157540.18) #dyrk1a
mean(normal$ENSG00000105204.12) #dyrk1b
tumoral <-merge(tumoral, rna, by="SAMPLE.ID")
mean(tumoral$ENSG00000105976.13) #c-met
mean(tumoral$ENSG00000157540.18) #dyrk1a
mean(tumoral$ENSG00000105204.12) #dyrk1b
t.test(tumoral$ENSG00000105976.13, normal$ENSG00000105976.13)[3] #c-met t-test 0.007631127
t.test(tumoral$ENSG00000157540.18, normal$ENSG00000157540.18)[3] #dyrk1a t-test ns
t.test(tumoral$ENSG00000105204.12, normal$ENSG00000105204.12)[3] #dyrk1b t-test ns
t.test(tumoral$ENSG00000146648.14, normal$ENSG00000146648.14)[3] #EGFR ns
t.test(tumoral$hsa.mir.148a, normal$hsa.mir.148a)[3] #mir-148a ns
tumoral2 <- tumoral[,c("DAYS.TO.DEATH", "AGE.AT.DIAGNOSIS", "ENSG00000157540.18")]
#tumoral2 <- tumoral2[complete.cases(tumoral2), ]
tumoral2 <- tumoral2[order(tumoral2$ENSG00000157540.18),] 
write.csv(x = tumoral2, "tumoral-dyrk1a.csv", row.names = FALSE)

mean(tumoral2[c(1:88),2])
sd(tumoral2[c(1:88),2])
mean(tumoral2[c(90:177),2])
sd(tumoral2[c(90:177),2])
t.test(tumoral2[c(1:88),2], tumoral2[c(90:177),2])[3]
t.test(tumoral2[c(1:44),2], tumoral2[c(134:177),2])[3]
