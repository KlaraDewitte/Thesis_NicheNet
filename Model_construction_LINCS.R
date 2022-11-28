library(nichenetr)
library(tidyverse)
library(dplyr)
# in the NicheNet framework, ligand-target links are predicted based on collected biological knowledge on ligand-receptor, signaling and gene regulatory interactions

# The complete networks can be downloaded from Zenodo
lr_network = readRDS(url("https://zenodo.org/record/5884439/files/lr_network_human_21122021.rds"))
sig_network = readRDS(url("https://zenodo.org/record/5884439/files/signaling_network_human_21122021.rds"))
gr_network = readRDS(url("https://zenodo.org/record/5884439/files/gr_network_human_21122021.rds"))
gr_network_LINCS = read.table("C:/Users/klara/MaStat/Thesis/epiNN/Ananse_networks/LINCS_MCF10A_TGFB_48_grn.txt", header = TRUE, sep = ";")

#add your new source to the weighted networks in NicheNet (weight??)
new_network_weights_df <- tibble(source = "LINCS_MCF10A_TGFB_48", weight = 1)
new_source_weights_df <- optimized_source_weights_df %>% bind_rows(new_network_weights_df)

#add LINCS GRN to default GRN NicheNet
new_gr_network <- gr_network %>% bind_rows(gr_network_LINCS)

# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = new_gr_network, source_weights_df = new_source_weights_df)

ligands <- list("EGF", "HGF", "OSM", "IFNG", "TGFB1", "BMP2")

ligand_target_matrix <- construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands)
