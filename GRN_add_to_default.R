library(nichenetr)
library(tidyverse)
library(dplyr)

#load background NicheNet networks
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))

#Path to where your network files are stored
path <- "C:/Users/klara/MaStat/Thesis/epiNN/Networks/"
files <- list.files(path = path, pattern = '.network', full.names = TRUE)

#add your source with weight to background NicheNet
new_network_weights_df <- tibble(source = "ReMap_2022", weight = 1)
new_source_weights_df <- optimized_source_weights_df %>% bind_rows(new_network_weights_df)

all_networks <- list()
all_names <- list()

for (file in files){
  name_file <- basename(file)
  all_names <- append(all_names, name_file)
  print(name_file)
  gr_network_remap <- read.csv(file, sep = '\t', col.names = c("from","pd","to"))
  gr_network_remap <- gr_network_remap %>% transmute(from = from,to = to) %>% mutate(source = "ReMap_2022", database = "ReMap")
  gr_network_remap <- gr_network_remap[-c(which(gr_network$to == "")),]
  new_gr_network <- gr_network %>% bind_rows(gr_network_remap)
  weighted_networks <- construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = new_gr_network, source_weights_df = new_source_weights_df)
  name <- strsplit(name_file, split = '_')
  ligands <- list(name[[1]][3])
  ligands <- gsub("IFNg", "IFNG", ligands)
  ligands <- gsub("TGFb", "TGFB1", ligands)
  ligands <- list(ligands)
  if(ligands != "PBS"){
    ligand_target_matrix <- construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands)
    weighted_networks$ligand_target_matrix <- as.data.frame(ligand_target_matrix)
    ligand_target_matrix <- tibble::rownames_to_column(as.data.frame(ligand_target_matrix), "VALUE")
    write.table(ligand_target_matrix, paste0("C:/Users/klara/MaStat/Thesis/epiNN/new_GRN", "/", name_file, "_ligand_target.matrix"), row.names = FALSE)
  }
  all_networks <- append(all_networks, list(weighted_networks))
}
names(all_networks) <- all_names


