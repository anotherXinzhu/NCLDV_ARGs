
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)


# NCLDV ======================================

# genome list -------------------
GenomeList <- read_excel("../_data/TableS.xlsx", sheet = 1, skip = 2)
GenomeList <- GenomeList[1:1416,]

# VF list ---------------------
VFlist <- read_excel("../_data/TableS.xlsx",sheet = 10, skip = 2)
VFlist[1,3:13]
colnames(VFlist)[3:13] <- as.character(VFlist[1,3:13])
VFlist <- VFlist[-1,]

all( VFlist$`Sequence ID` %in% GenomeList$`Sequence ID a` )

GenomeList$VF <- 
  sapply(GenomeList$`Sequence ID a`,
         function(x) if(x %in% VFlist$`Sequence ID`) "VF+" else "VF-")

plotDat_ncldv <- 
  table(GenomeList$VF) %>%
  as.data.frame() %>%
  mutate(totalN = sum(Freq)) %>%
  mutate(perc = Freq/totalN)

ggplot(plotDat_ncldv) +
  geom_col(aes(x=Var1, y=perc, fill=Var1), width = 0.6, color="black") +
  scale_fill_manual(values = c("#FFF8BC","#6897C2")) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  ylab("Genomes (%)") + xlab("")
ggsave("Fig6A.VF_GenomePerc_NCLDV.pdf", width = 3, height = 2.2)




# phage ================================================


# genome list -------------------
GenomeList <- read_excel("../_data/TableS.xlsx", sheet = 2,   skip = 2)
GenomeList <- GenomeList[1:40277,]

# VF list ---------------------
VFlist <- read_excel("../_data/TableS.xlsx", sheet = 11, skip = 2)
VFlist[1,3:13]
colnames(VFlist)[3:13] <- as.character(VFlist[1,3:13])
VFlist <- VFlist[-1,]


all( VFlist$`Sequence ID` %in% GenomeList$`Sequence ID a` )

GenomeList$VF <- 
  sapply(GenomeList$`Sequence ID a`,
         function(x) if(x %in% VFlist$`Sequence ID`) "VF+" else "VF-")

plotDat_phage <- 
  table(GenomeList$VF) %>%
  as.data.frame() %>%
  mutate(totalN = sum(Freq)) %>%
  mutate(perc = Freq/totalN)


ggplot(plotDat_phage) +
  geom_col(aes(x=Var1, y=perc, fill=Var1), width = 0.6, color="black") +
  scale_fill_manual(values = c("#FFF8BC","#6897C2")) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  ylab("Genomes (%)") + xlab("")
ggsave("Fig6E.VF_GenomePerc_phage.pdf",device = "pdf", width = 3, height = 2.2)
