library(phyloseq)
library(magrittr)
setwd('F:/DataPaula/Esparanza/')


taxa <-  read.csv2("Taxonomy.txt", sep = '\t', row.names = 1, stringsAsFactors = F) %>% as.matrix
taxa[which(taxa=="")] <- NA

taxa %<>% tax_table()

asv <- read.csv2("filtered_table.txt", sep = '\t', row.names = 1, header = T, stringsAsFactors = F) %>% otu_table(taxa_are_rows = T)

meta <- read.csv2("Metadata.csv", stringsAsFactors = F, row.names = 1) %>% sample_data()


asv_table <- phyloseq(taxa, asv, meta)

asv_table.fam <- tax_glom(asv_table, taxrank = 'Family')
asv_table.fam.relab <- transform_sample_counts(asv_table.fam, function(OTU) OTU/sum(OTU))