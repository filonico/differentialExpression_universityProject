library(NOISeq)
library(tidyverse)

#upload raw counts
dgal <- data.frame(read.table("all_rawcounts_filtered.txt", header = T, row.names = 1))
conditions <- data.frame(dgal_pool = c("M6_C", "M6_C", "M6_C", "M6_F", "M6_F", "M6_F"))

#eSet objet  
dgal_eset <- readData(data = dgal, factors = conditions)

###DATA QUALITY ANALYSIS###
#saturation plot
dgal_eset %>% dat(k = 0, ndepth = 7, type = "saturation") %>% explo.plot(toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)

#raw count barplot
dgal_eset %>% dat(factor = NULL, type = "countsbio") %>% explo.plot(toplot = 1, samples = NULL, plottype = "barplot")

#filtering reads and relative eSet object
dgal_filtered <- filtered.data(dgal, factor = conditions$dgal_pool, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 10, p.adj = "fdr")
dgal_filtered_eset <- readData(data = dgal_filtered, factors = conditions)

#filtered saturation plot
dgal_filtered_eset %>% dat(k = 0, ndepth = 7, type = "saturation") %>% explo.plot(toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)

#filtered raw count barplot
dgal_filtered_eset %>% dat(factor = NULL, type = "countsbio") %>% explo.plot(toplot = 1, samples = NULL, plottype = "barplot")

#write normalization table
assayData(dgal_filtered_eset)$exprs %>% tmm(long = 1000, lc = 0) %>% write.table(file="dgal_tmm10.txt", append = FALSE, eol="\n", quote = FALSE)

#differential expression
assayData(dgal_filtered_eset)$exprs %>% tmm(long = 1000, lc = 0) %>% readData(factors = conditions) %>% noiseqbio(k=0.1, factor="dgal_pool") %>% degenes(q = 0.95, M = NULL) %>% write.table(file="dgal_diff_expr.txt", append = FALSE, eol="\n", quote = FALSE)

#plot the average expression value and highlight the feature differentially expressed
assayData(dgal_filtered_eset)$exprs %>% tmm(long = 1000, lc = 0) %>% readData(factors = conditions) %>% noiseqbio(k=0.1, factor="dgal_pool") %>% DE.plot(q = 0.95, graphic = "expr", log.scale = TRUE)

assayData(dgal_filtered_eset)$exprs %>% tmm(long = 1000, lc = 0) %>% readData(factors = conditions) %>% noiseqbio(k=0.1, factor="dgal_pool") %>% DE.plot(q = 0.95, graphic = "MD")

#plot differential expression
diff_expr <- data.frame(read.table("dgal_diff_expr.txt", header=T, dec="."))

#diff_expr[,c(1,2,3)] %>% pivot_longer(-1) %>% ggplot(aes(geneID, value, fill=name)) + 
# geom_bar(position="dodge", stat="identity") +
#  scale_y_log10()
