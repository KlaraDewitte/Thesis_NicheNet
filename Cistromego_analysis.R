library(tidyverse)
library(stringr)

path <- "C:/Users/klara/MaStat/Thesis/GithubThesis/Thesis_21_22/Files"
files <- list.files(path = path, pattern = '.csv', full.names = TRUE)

csv_files <- list()
all_names <- list()

for (file in files){
  name <- basename(file)
  values <- stringr::str_split(name, "_") %>% 
            do.call(rbind, .) 
  unique_values <- values %>%
                   apply(2, unique)          
  val <- paste0(unique_values[3], "_", unique_values[4], "_", unique_values[6])
  k <- read.csv(file, header = FALSE, sep='\t') %>% 
      rename(Gene = V1, Score = V2)
  rank <- rank(-k[,2], ties.method = "min")
   
  l <- cbind(k, rank)
  l <- list(l[order(l$Gene),])
  all_names <- append(all_names, val)
  csv_files <- append(csv_files, l)
}
names(csv_files) <- all_names

values <- stringr::str_split(all_names, "_") %>% 
            do.call(rbind, .)
unique_values <- values %>%
                  apply(2, unique) %>%
                  setNames(c("cytokine", "timepoint", "replicate"))

all_names_2 <- list()
for(ctx in unique_values$cytokine){
  for(tmp in unique_values$timepoint){
    all_names_2 <- append(all_names_2, paste0(ctx, "_", tmp))
    if(any((str_detect(all_names, paste0(ctx, "_", tmp))==TRUE))){
      genes <- as.data.frame(csv_files[[1]][["Gene"]])
      total_rank <- data.frame(genes)
      total_scores <- data.frame(genes)
      colnames(total_rank) <- 'Gene'
      colnames(total_scores) <- 'Gene'
      number <- which(str_detect(all_names, paste0(ctx, "_", tmp), negate = FALSE))
      for(l in number){
        total_rank <- cbind(total_rank, csv_files[[l]][["rank"]])
        total_scores <- cbind(total_scores, csv_files[[l]][["Score"]])
        names(total_rank)[length(names(total_rank))] <- all_names[l]
        names(total_scores)[length(names(total_scores))] <- all_names[l]
      }
      median_rank <- apply(total_rank[,-1], 1, median)
      ranks <- order(median_rank)
      total_rank$row_median <-  ranks
      median_score <- apply(total_scores[,-1], 1, median)
      total_scores$row_median_score <- median_score
      score_ranks <- rank(-total_scores[,"row_median_score"], ties.method = "min")
      score <- data.frame(genes)
      colnames(score) <- 'gene'
      score$score <- median_score
      total_scores$row_rank <- score_ranks
      score$rank <- score_ranks
      write.csv(total_rank, paste0("C:/Users/klara/MaStat/Thesis/GithubThesis/Thesis_21_22/Files/Analysis", "/", paste0(ctx, "_", tmp), "_rank.csv"), row.names = FALSE)
      write.csv(total_scores, paste0("C:/Users/klara/MaStat/Thesis/GithubThesis/Thesis_21_22/Files/Analysis", "/", paste0(ctx, "_", tmp), "_scores.csv"), row.names = FALSE)
      write.csv(score, paste0("C:/Users/klara/MaStat/Thesis/GithubThesis/Thesis_21_22/Files/Analysis/Analysis", "/", paste0(ctx, "_", tmp), ".csv"), row.names = FALSE)
       }
    }
}





