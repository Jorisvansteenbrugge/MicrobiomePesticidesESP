library(phyloseq)
library(magrittr)
setwd('F:/DataPaula/Esparanza/')

## export otu ids for phylogeny
#asv_table.fam %>% taxa_names() %>% cat

#import tree for phylogenic based analyses
tree <- ape::read.tree("fasta_phylogeny/otu_sequences_subset_MSA.tree") %>% phy_tree()




taxa <-  read.csv2("Taxonomy.txt", 
                   sep = '\t', 
                   row.names = 1, 
                   stringsAsFactors = F) %>% as.matrix

taxa[which(taxa=="")] <- NA

taxa %<>% tax_table()

asv <- read.csv2("filtered_table.txt", 
                 sep = '\t', 
                 row.names = 1, 
                 header = T, 
                 stringsAsFactors = F) %>% otu_table(taxa_are_rows = T)

meta <- read.csv2("Metadata.csv", 
                  stringsAsFactors = F, 
                  row.names = 1) %>% sample_data()

meta_treatment <- read.csv2("Meta_treatment.csv", 
                            stringsAsFactors = F, 
                            row.names = 1) %>% sample_data()
  

  

asv_table <- phyloseq(taxa, 
                      asv, 
                      meta)

asv_table.fam <- tax_glom(asv_table, 
                          taxrank = 'Family')
phy_tree(asv_table.fam) <- tree

asv_table.fam.relab <- transform_sample_counts(asv_table.fam, function(OTU) OTU/sum(OTU))


asv_table.treatment <- phyloseq(taxa, 
                                asv, 
                                meta_treatment)

asv_table.treatment.fam <- tax_glom(asv_table.treatment, 
                                    taxrank = 'Family')
phy_tree(asv_table.treatment.fam) <- tree
asv_table.treatment.fam.relab <- transform_sample_counts(asv_table.treatment.fam, function(OTU) OTU/sum(OTU))


save(asv_table.fam, 
     asv_table.fam.relab, 
     asv_table.treatment.fam, 
     asv_table.treatment.fam.relab, 
     file = "Asv_tables.RData")





