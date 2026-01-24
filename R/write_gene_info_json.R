library(jsonlite)
library(here)
library(tmod)
library(dplyr)
library(tidyr)
library(magrittr)
library(readr)
library(ggplot2)
library(readxl)
library(tidyr)
#library(GO.db)
#library(org.Hs.eg.db)
library(here)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")

#data from:
#curl -L "https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606)+AND+(reviewed:true)&format=tsv&fields=accession,id,gene_primary,protein_name,cc_function" \
#> uniprot_swissprot_human_function.tsv
gene_data <- read_tsv(here("from_uniprot/uniprot_swissprot_human_function.tsv"))

target_cluster_path <- here("v6")
cluster_mapping <- read_csv(here(target_cluster_path, "HBSCAN_clusters.csv")) 

gene_data <- left_join(cluster_mapping %>% select(gene_symbol), gene_data %>% select(gene_symbol = `Gene Names (primary)`, uniprot_entry=Entry, name = `Protein names`, description =  `Function [CC]` ))

gene_data %<>% mutate(description = gsub("FUNCTION: ", "", description))



gene_json <- gene_data %>%
  transmute(
    gene_symbol,
    entry = Map(function(n, d, u) {
      base <- list(name = n, uniprot = u)
      if (is.na(d) || d == "") base else c(base, list(description = d))
    }, name, description, uniprot_entry)
  ) %>%
  deframe()

write_json(gene_json, here(target_cluster_path, "gene-metadata.json"), auto_unbox = TRUE, pretty = TRUE)
