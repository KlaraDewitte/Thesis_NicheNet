library(nichenetr)
library(tidyverse)
library(dplyr)

#load background NicheNet networks
lr_network = readRDS(url("https://zenodo.org/record/5884439/files/lr_network_human_21122021.rds"))
sig_network = readRDS(url("https://zenodo.org/record/5884439/files/signaling_network_human_21122021.rds"))

#Path to where your network files are stored
path <- "C:/Users/klara/MaStat/Thesis/epiNN/Networks/"
files <- list.files(path = path, pattern = '.network', full.names = TRUE)

#add your source with weight to background NicheNet
source_weights_df <- add_row(source_weights_df, source = "ReMap_2022", weight = 1)

all_networks <- list()
all_names <- list()

for (file in files){
  name_file <- basename(file)
  all_names <- append(all_names, name_file)
  print(name_file)
  gr_network <- read.csv(file, sep = '\t', col.names = c("from","pd","to"))
  gr_network <- subset(gr_network, select = -c(pd)) %>% cbind(rep('ReMap_2022', nrow(gr_network)))
  colnames(gr_network) <- c("from", "to", "source")
  gr_network <- gr_network[-c(which(gr_network$to == "")),]
  weighted_networks <- construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, source_weights_df = source_weights_df)
  name <- strsplit(name_file, split = '_')
  ligands <- list(name[[1]][3])
  ligands <- gsub("IFNg", "IFNG", ligands)
  ligands <- gsub("TGFb", "TGFB1", ligands)
  ligands <- list(ligands)
  if(ligands != "PBS"){
    ligand_target_matrix <- construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands)
    weighted_networks$ligand_target_matrix <- as.data.frame(ligand_target_matrix)
    ligand_target_matrix <- tibble::rownames_to_column(as.data.frame(ligand_target_matrix), "VALUE")
    write.table(ligand_target_matrix, paste0("C:/Users/klara/MaStat/Thesis/epiNN/ligand_target_matrix", "/", name_file, "_ligand_target.matrix"), row.names = FALSE)
  }
  all_networks <- append(all_networks, list(weighted_networks))
}
names(all_networks) <- all_names


