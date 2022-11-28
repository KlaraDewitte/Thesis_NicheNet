library(nichenetr)
library(tidyverse)
library(dplyr)
# in the NicheNet framework, ligand-target links are predicted based on collected biological knowledge on ligand-receptor, signaling and gene regulatory interactions

# The complete networks can be downloaded from Zenodo
lr_network = readRDS(url("https://zenodo.org/record/5884439/files/lr_network_human_21122021.rds"))
sig_network = readRDS(url("https://zenodo.org/record/5884439/files/signaling_network_human_21122021.rds"))
gr_network = readRDS(url("https://zenodo.org/record/5884439/files/gr_network_human_21122021.rds"))

matrix <- read.csv("C:/Users/klara/MaStat/Thesis/epiNN/Matrices/LIB180727LH_02_TGFb_24_S48_L006__001_ATAC_TF_TG.matrix", sep= "\t")

test <- read.csv("C:/Users/klara/MaStat/Thesis/epiNN/Networks/LIB180727LH_02_TGFb_24_S48_L006__001_ATAC.network", sep ="\t", col.names = c("from","pd","to"))
test <- subset(test, select = -c(pd)) %>% cbind(rep('ReMap_2022', nrow(test))) 
colnames(test) <- c("from", "to", "source")
testtest <- test[-c(which(test$to == "")),]
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = testtest, source_weights_df = source_weights_df)

ligand_target_matrix <- construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands)

colnames(test) <- c("from", "to", "source")

path <- "C:/Users/klara/MaStat/Thesis/epiNN/Matrices/"
files <- list.files(path = path, pattern = '.matrix', full.names = TRUE)

for (file in files){
  name <- basename(file)
  print(name)
  matrix <- read.csv(file, sep = '\t')
  gr <- data.frame()
  for(r in 1:nrow(matrix)){
    o <- rep(matrix$GENE[r], length(colnames(matrix)[which(matrix[r,]==1)]))
    network <- data.frame(o)
    network <- cbind(network, colnames(matrix)[which(matrix[r,]==1)])
    gr <- rbind(gr, network)
  }
  gr <- cbind(gr, rep('ReMap_2022', nrow(gr)))
  colnames(gr) <- c("to", "from", "source")
  write.table(gr, paste0("C:/Users/klara/MaStat/Thesis/epiNN/Matrices/Analysis", "/", name, "_gr.matrix"), row.names = FALSE)
}

gr <- read.table("C:/Users/klara/MaStat/Thesis/epiNN/Matrices/Analysis/LIB180727LH_02_HGF_24_S24_L004__001_ATAC_TF_TG.matrix_gr.matrix", header = TRUE)
gr_new <- gr[-c(which(gr$from == "")),]
colnames(gr_new) <- c("to", "from", "source")

lr_network = readRDS(url("https://zenodo.org/record/5884439/files/lr_network_human_21122021.rds"))
sig_network = readRDS(url("https://zenodo.org/record/5884439/files/signaling_network_human_21122021.rds"))
path <- "C:/Users/klara/MaStat/Thesis/GithubThesis/Thesis_21_22/Matrices/Analysis"
files <- list.files(path = path, pattern = '.matrix', full.names = TRUE)
source_weights_df <- add_row(source_weights_df, source = "ReMap_2022", weight = 1)

all_weighted_networks <- list()
all_names <- list()
for (file in files){
  name <- basename(file)
  all_names <- append(all_names, name)
  print(name)
  gr_network <- read.table(file, header = TRUE)
  gr_network <- gr_network[-c(which(gr_network$from == "")),]
  weighted_networks <- construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = test, source_weights_df = source_weights_df)
  all_weighted_networks <- append(all_weighted_networks, list(weighted_networks))
}
names(all_weighted_networks) <- all_names
for(i in 1:length(all_weighted_networks)){
  weighted_networks <- all_weighted_networks[[i]]
  name <- names(all_weighted_networks[i])
  name <- strsplit(name, split = '_')
  ligands <- list(name[[1]][3])
  if (ligands == "TGFB"){
    ligands <- list("TGFB1")
  }
  if(ligands != "PBS"){
    ligand_target_matrix <- construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands)
    all_weighted_networks[[i]]$ligand_target_matrix <- ligand_target_matrix
  }
}


for (i in 1:length(all_networks)){
  all_networks[[i]]$ligand_target_matrix <- as.data.frame(all_networks[[i]]$ligand_target_matrix)
  all_networks[[i]]$ligand_target_matrix <- tibble::rownames_to_column(all_networks[[i]]$ligand_target_matrix, "VALUE")
}



ligands <-  list(unique_values$cytokine[-1:-2])
ligands <- ligands[[1]][-3]
ligands <- append(ligands, "TGFB1")
ligands <- as.list(ligands)

ligands <- list("HGF")
ligands[[1]]
ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands)

