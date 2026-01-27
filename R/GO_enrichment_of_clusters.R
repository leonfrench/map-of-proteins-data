library(future)
library(future.apply)
library(tmod)
library(dplyr)
library(tidyr)
library(magrittr)
library(readr)
library(ggplot2)
library(readxl)
library(tidyr)
library(GO.db)
library(org.Hs.eg.db)
library(here)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")

target_cluster_path <- here("v6")
cluster_mapping <- read_csv(here(target_cluster_path, "HBSCAN_clusters.csv"))

go_table <- as.data.frame(org.Hs.egGO2ALLEGS) #gene mapping to the direct GO term and that term's children
go_table %<>% as_tibble()
go_table %<>% mutate(gene_symbol = mapIds(org.Hs.eg.db, gene_id, "SYMBOL","ENTREZID")) #Add gene symbol
go_table <- go_table %>% select(gene_symbol, go_id) #we only need the sybmol and the GO ID at this stage
go_table %<>% distinct()
go_table %<>% filter(gene_symbol %in% cluster_mapping$gene_symbol)
go_table %<>% mutate(GO_name = Term(go_id))

#create tmod object
tmodNames <- tibble()
highlevel_map <- list()
high_level_size_threshold <- 100
low_level_size_threshold <- 5
right_sized_GO <- go_table %>% group_by(GO_name) %>% count() %>% filter(n < high_level_size_threshold, n> low_level_size_threshold) %>% pull(GO_name)

#create a tmod object
for(high_level_name in right_sized_GO) {
  if(nrow(tmodNames) %% 1000 == 0) { print(nrow(tmodNames))}
  #print(high_level_name)
  gene_symbols <- unique(go_table %>% filter(GO_name == high_level_name) %>% pull(gene_symbol))  
  if(length(gene_symbols) > high_level_size_threshold) {  #filter for greater than a threshold
    next
  }
  if(length(gene_symbols) < low_level_size_threshold) {  #filter for greater than a threshold
    next
  }
  highlevel_map[high_level_name] <- list(gene_symbols)
  tmodNames <- bind_rows(tmodNames, data.frame(ID=high_level_name, Title = high_level_name))
}
tmodNames
GO_tmod_obj <- makeTmod(modules = tmodNames, modules2genes = highlevel_map)
GO_tmod_obj


background_genes <- cluster_mapping$gene_symbol

# ---- parallel HG test per cluster ----
background_genes <- cluster_mapping$gene_symbol
go_genes <- unique(go_table$gene_symbol)

# Pre-split once so we don't scan cluster_mapping over and over
cluster_genes_map <- split(cluster_mapping$gene_symbol, cluster_mapping$hdbscan_filled_cluster)
cluster_genes_map <- lapply(cluster_genes_map, function(gs) gs[gs %in% go_genes])

clusters <- names(cluster_genes_map)

# Optional restart: if an output already exists, skip finished clusters
out_path <- here(target_cluster_path, "gene_cluster_xy.GO_enrich.csv")

# Use most cores minus 1 (tweak if you want)
workers <- max(1, future::availableCores() - 2)
future::plan(future::multisession, workers = workers)

res_list <- future.apply::future_lapply(clusters, function(target_cluster) {
  genes_in_cluster <- cluster_genes_map[[target_cluster]]
  
  result <- tmod::tmodHGtest(
    fg = genes_in_cluster,
    bg = background_genes,
    mset = GO_tmod_obj,
    qval = 1.01,
    filter = FALSE
  )
  result <- tibble::as_tibble(result)
  
  if (nrow(result) == 0) {
    return(tibble::tibble(hdbscan_filled_cluster = target_cluster, GO_hit = "No significant GO enrichment"))
  }
  
  result <- dplyr::filter(result, adj.P.Val < 0.1)
  if (nrow(result) == 0) {
    return(tibble::tibble(hdbscan_filled_cluster = target_cluster, GO_hit = "No significant GO enrichment"))
  }
  
  tibble::tibble(hdbscan_filled_cluster = target_cluster, GO_hit = result$ID[1])
})
#should just be a matrix multiplication
future::plan(future::sequential)

all_HG_results <- dplyr::bind_rows(res_list)

all_HG_results <- all_HG_results %>% dplyr::mutate(hdbscan_filled_cluster = as.integer(hdbscan_filled_cluster))

all_HG_results %>% pull(GO_hit) %>% n_distinct()
all_HG_results %>% group_by(GO_hit) %>% count() %>% arrange(-n) %>% filter(n>1)
all_HG_results %<>% inner_join(cluster_mapping %>% group_by(hdbscan_filled_cluster) %>% summarize(gene_symbols = paste(gene_symbol, sep = " ", collapse = " ")))
all_HG_results %<>% mutate(annotated_label = "")

all_HG_results %<>% arrange(GO_hit)
all_HG_results %>% write_csv(here(target_cluster_path, "gene_cluster_xy.GO_enrich.csv")) 
here(target_cluster_path, "gene_cluster_xy.GO_enrich.csv")


#ChatGPT prompt that was used on deep research with 1/3 of the clusters each time (ID, GO group and genes as input):
# You are a master namer of protein clusters. Your task is to create a unique, specific, and memorable name for a set of human protein coding genes based on their common function or subcellular localization.
# 
# The name should be:
# 
# 1. Concise (1-3 words maximum)
# 2. Distinctive and specific to these particular proteins/genes
# 3. Capture the unique essence of this specific set
# 4. AVOID generic terms
# 5. Creative but immediately understandable
# 6. If proteins are focused on a specific function or location, the name should reflect this specificity
# 7. IMPORTANT: If two protein cluster names are similar, DO NOT combine their names (e.g., avoid "NodeNexus" if there's also a "Node Nexus")
# 8. Use strong imagery and metaphors
# Each name must be HIGHLY DISTINCT from all other protein cluster names. Imagine this name appearing on a map - it should be instantly recognizable.
# 
# Don't write code, inspect it manually. Double-check your work.
# Do this for all of the clusters in the below table. This table gives the most enriched Gene Ontology group and the names of the genes for reference. The Gene Ontology group name can guide your name but don't copy it completely. For output just give the cluster ID and the name.



