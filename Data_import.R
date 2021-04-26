library(phyloseq)
library(magrittr)
setwd('F:/DataPaula/Esparanza/')


taxa <-  read.csv2("Taxonomy.txt", sep = '\t', row.names = 1, stringsAsFactors = F) %>% as.matrix
taxa[which(taxa=="")] <- NA

taxa %<>% tax_table()

asv <- read.csv2("filtered_table.txt", sep = '\t', row.names = 1, header = T, stringsAsFactors = F) %>% otu_table(taxa_are_rows = T)

meta <- read.csv2("Metadata.csv", stringsAsFactors = F, row.names = 1) %>% head %>% sample_data()


asv_table <- phyloseq(taxa, asv, meta)
