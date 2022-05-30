if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("NOISeq")

install.packages("gridExtra")

library(tidyverse)
library(gridExtra)

rc_mapped = tibble(read.table("all_rawcounts_mapped.tsv", header=T))

plot_mapped <- rc_mapped %>% gather(key="conditions", value="rawcounts", -1) %>% mutate(log_rawcounts = log(rawcounts)) %>%
  ggplot(aes(geneID, conditions, fill = log_rawcounts)) + 
  geom_tile() +
  scale_fill_distiller(palette = "BuGn", direction = +1) +
  labs(y = "Runs", fill = "#reads\n(log)") + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

rc_unmapped = tibble(read.table("all_rawcounts_unmapped.tsv", header=T))

plot_unmapped <- rc_unmapped %>% gather(key="conditions", value="rawcounts", -1) %>% mutate(log_rawcounts = log(rawcounts)) %>%
  ggplot(aes(geneID, conditions, fill = log_rawcounts)) + 
  geom_tile() +
  scale_fill_distiller(palette = "OrRd", direction = +1) +
  labs(x = "Transcripts", y = "Runs", fill = "#reads\n(log)") + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

grid.arrange(plot_mapped, plot_unmapped, nrow = 2)

rc_filtered = tibble(read.table("all_rawcounts_filtered.txt", header=T))

plot_filtered <- rc_filtered %>% gather(key="conditions", value="rawcounts", -1) %>% mutate(log_rawcounts = log(rawcounts)) %>%
  ggplot(aes(geneID, conditions, fill = log_rawcounts)) + 
  geom_tile() +
  scale_fill_distiller(palette = "BuGn", direction = +1) +
  labs(x = "Transcripts", y = "Runs", fill = "#reads\n(log)") + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
