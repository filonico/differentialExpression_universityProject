library(topGO)
library(tidyverse)

#create topGO file of GO annotation
dgal_ID2GO <- readMappings(file = "pannzer_annotation_rd_topgo_tot.tsv", sep = "\t", IDsep = ",")

#get genes of interest
dgal_geneinterest <- read.table("dgal_diff_expr.txt", header=T, dec=".")

#get genes and the relative q-values
dgal_geneinterest_pval <- dgal_geneinterest[,c(1,5)]

#get FDRs from q-values
dgal_geneinterest_pval$prob <- as.vector(1-dgal_geneinterest_pval$prob)

#convert the gene and FDR dataframe into a named vector
dgal_geneinterest_pval_vector <- dgal_geneinterest_pval$prob
names(dgal_geneinterest_pval_vector) <- dgal_geneinterest_pval$geneID

#declare function to select the most differentially expressed genes
topDiffGenes <- function(allScore) {return(allScore < 0.01)}

#create the topGOdata object for each GO major class (namely, MF, BP and CC)
dgal_GOdataBP <- new("topGOdata", ontology = "BP", allGenes = dgal_geneinterest_pval_vector, geneSelectionFun = topDiffGenes, annot = annFUN.gene2GO, gene2GO = dgal_ID2GO)

dgal_GOdataCC <- new("topGOdata", ontology = "CC", allGenes = dgal_geneinterest_pval_vector, geneSelectionFun = topDiffGenes, annot = annFUN.gene2GO, gene2GO = dgal_ID2GO)

dgal_GOdataMF <- new("topGOdata", ontology = "MF", allGenes = dgal_geneinterest_pval_vector, geneSelectionFun = topDiffGenes, annot = annFUN.gene2GO, gene2GO = dgal_ID2GO)

#perform GO enrichment tests
dgal_kstest_BP <- runTest(dgal_GOdataBP, algorithm = "classic", statistic = "ks")

dgal_kstest_CC <- runTest(dgal_GOdataCC, algorithm = "classic", statistic = "ks")

dgal_kstest_MF <- runTest(dgal_GOdataMF, algorithm = "classic", statistic = "ks")

#write results into a tsv file 
dgal_kstest_BP %>% GenTable(classicKS = dgal_kstest_BP, ranksOf = "classic", topNodes = 35) %>% write.table(file="dgal_topGO_BP.tsv", quote=FALSE, row.names=FALSE, sep = "\t")

dgal_kstest_CC %>% GenTable(classicKS = dgal_kstest_CC, ranksOf = "classic", topNodes = 10) %>% write.table(file="dgal_topGO_CC.tsv", quote=FALSE, row.names=FALSE, sep = "\t")

dgal_kstest_MF %>% GenTable(classicKS = dgal_kstest_MF, ranksOf = "classic", topNodes = 18) %>% write.table(file="dgal_topGO_MF.tsv", quote=FALSE, row.names=FALSE, sep = "\t")
