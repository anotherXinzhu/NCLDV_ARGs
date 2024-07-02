
library(dplyr)
library(readxl)

dat1 <- read_excel("GX_tailings_Scaff0041861.annotation.xlsx", sheet = 1) %>% select(-1,-2)
dat2 <- read_excel("GCA_002116175.1.annotation.xlsx", sheet = 1 ) %>% select(-1,-2)

plotDat <- bind_rows(dat1 %>% select(scaffold, gene_position, start_position, end_position, strandedness, Classification, Text),
                     dat2 %>% select(scaffold, gene_position, start_position, end_position, strandedness, Classification, Text))
Classes <- plotDat$Classification %>% unique

library(paletteer)
paletteer_d("DresdenColor::colddays")

tmp1 <- cbind.data.frame(
  Class = Classes[grep("[Mm]etabolism", Classes)],
  Color = c('#D3EEF9FF',"#ABC9C8FF", "#4692B0FF", "#5074A5FF", "#2F70A1FF")
)

Classes[!grepl("[Mm]etabolism", Classes)]
tmp2 <- cbind.data.frame(
  Class = c("Translation","Transcription",  "DNA modification or repair","Drug resistance","Others","Unidentified","Unknown function"),
  Color = c("#7F4B89FF", "#DFE07CFF", "#E3A5D6FF","#C00000FF","#BF714DFF","white", "white")
)

Colors_df <- bind_rows(tmp1, tmp2)
Colors <- Colors_df$Color
names(Colors) <- Colors_df$Class

plotDat$orientation = sapply(plotDat$strandedness, function(x) if(x == 1) 1 else 0)

library(ggplot2)
library(gggenomes)
pdf("geneStructure.pdf", width = 12,height = 3)
ggplot(plotDat,
       aes(xmin = start_position, xmax = end_position, y = scaffold, fill = Classification, 
           forward = orientation)) +
  geom_gene_arrow() +
  facet_wrap(~ scaffold, scales = "free", ncol = 1) +
  theme_genes() +
  scale_fill_manual(values = Colors) +
  theme(legend.position="bottom") +
  geom_text(aes(x=(start_position + end_position)/2, y=scaffold , label=Text, vjust = -1))
dev.off()


