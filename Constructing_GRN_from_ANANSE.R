library(dplyr)
library(stringr)

timepoints <- c(24, 48)
ligands <- c("TGFB", "OSM", "EGF", "HGF")


for(t in timepoints){
  for(l in ligands){
    print(paste0("Ananse_networks/LINCS_MCF10A_", l, "_", t, "_NO_RNAseq.network.txt"))
    line1 <- read.delim(file = paste0("Ananse_networks/LINCS_MCF10A_", l, "_", t, "_NO_RNAseq.network.txt"), header = T) %>%
      filter(prob >= 0.8) %>% mutate_all(~gsub("Râ€”", "-", .)) %>% mutate_all(~gsub("â€”", "-", .)) %>%
      subset(, select = -c(prob))
    
    print((paste0("LINCS_MCF10A_", l, "_", t)))
    
    tf_target_column <- subset(line1, select = c("tf_target"))
    
    for(k in 1:nrow(tf_target_column)){
      splitted <- str_split(tf_target_column[k,], "-", n = Inf, simplify = FALSE)
      line1[k, "from"] <- splitted[[1]][1]
      line1[k, "to"] <- splitted[[1]][2]
      print(k)
    }
    
    line1 <- subset(line1, select = -c(tf_target))  %>%
      cbind(paste0("LINCS_MCF10A_", l, "_", t)) %>% cbind("LINCS_MCF10A") 
    colnames(line1) <- c("to", "from", "source", "database")
    
    write.csv(line1, paste0("Ananse_networks/GRN/LINCS_MCF10A_", l, "_", t, ".csv"), row.names=FALSE)
    gc()
  }
}


