
library(dplyr)
library(readxl)
library(Biostrings)
TableS2 <- read_excel("2019-Murina-TableS2_all_ABCF_sequences.xlsx")
taxonomy <- read_excel("2019-Murina-TableS1_Taxonomy_ABCF_composition.xlsx")

TableS2$Superkingdom <- sapply(TableS2$Species, function(x) taxonomy$Superkingdom[which(taxonomy$Species == x)])  
TableS2$seqLabel = paste(TableS2$Superkingdom,TableS2$Subfamily, seq(1,nrow(TableS2),1), sep = "_")

library(Biostrings)
test <- AAStringSet(TableS2$Sequence) #turn characterinto biostring object
names(test) <- TableS2$seqLabel
writeXStringSet(test, filepath = 'ABCF_refSeqs.faa')
write.table(TableS2 %>% select(-Sequence), 
            file = "meta.ABCF_refSeqs.txt",
            sep = "\t", quote = F, row.names = F)

# randomly select 5 sequences for each subfamily
set.seed(1)
test <- TableS2 %>%
  filter(Subfamily != "ABCF1/2") %>%  # ABCF1/2 only has two sequences
  dplyr::group_by(Subfamily) %>%
  dplyr::summarise(seqLabel = sample(seqLabel, size=5)) %>%
  as.data.frame()
refSeq.Selected <- bind_rows(test, TableS2 %>% filter(Subfamily == "ABCF1/2") %>% select(Subfamily, seqLabel)) 
sapply(refSeq.Selected, class)


# combine sequences to construct tree
refSeqs <- readAAStringSet("ABCF_refSeqs.faa")
viralSeqs <- readAAStringSet("ABCF_viralSeqs.faa")
treeSeqs <- c(refSeqs[refSeq.Selected$seqLabel], viralSeqs)

writeXStringSet(treeSeqs, filepath = "ABCF_treeSeqs.faa")