
library(data.table)
library(readxl)


# NCLDV =====================================================
# genome list -------------------
GenomeList <- read_excel("../_data/TableS.xlsx", sheet = 1, skip = 2)
GenomeList <- GenomeList[1:1416,]


# VF list ---------------------
VFlist <- read_excel("../_data/TableS.xlsx",sheet = 10, skip = 2)
VFlist[1,3:13]
colnames(VFlist)[3:13] <- as.character(VFlist[1,3:13])
VFlist <- VFlist[-1,]

# ARG list ----------------------
ARGlist <- read_excel("../_data/TableS.xlsx", sheet=4, skip = 2)
as.character(ARGlist[1,4:6])
colnames(ARGlist)[4:6] <- as.character(ARGlist[1,4:6])
ARGlist <- ARGlist[2:750,]


# 1）chisq.test  ---------
tmp <- 
  merge(
    GenomeList %>% select(`Sequence ID a`),
    ARGlist %>% mutate(ARG = "ARG+") %>% select(`Sequence ID`,ARG),
    by.x = "Sequence ID a", by.y = "Sequence ID",
    all.x = T)
tmp$ARG[is.na(tmp$ARG)] <- "ARG-"

tmp <- 
  merge(
    tmp,
    VFlist %>% mutate(VF="VF+") %>% select(`Sequence ID`, VF),
    by.x = 'Sequence ID a', by.y = "Sequence ID", all.x = T
  )
tmp$VF[is.na(tmp$VF)] <- "VF-"

table(tmp$ARG, tmp$VF)
ct <- chisq.test(tmp$ARG, tmp$VF)
ct$statistic
ct$p.value 


# 2） bar plot  ---------
library(ggplot2)

plotDat_ncldv <- table(tmp$VF, tmp$ARG) %>%
  as.data.frame() %>%
  group_by(Var1) %>%
  mutate(bar.total = sum(Freq)) %>%
  as.data.frame() %>%
  mutate(perc = Freq/bar.total)


ggplot(plotDat_ncldv %>% filter(Var2 == "ARG+") ) +
  geom_col(aes(x=Var1, y=perc,fill = Var1),width = 0.6, color="black") +
  theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#FFF8BC","#6897C2"))  +
  xlab("") + ylab("ARG-carrying genomes (%)")

ggsave(filename = "Fig6B.VF-ARG_chiSquare_ncldv.pdf", width = 3, height = 2.2, device = "pdf")





# phage ====================================================

# genome, ARG and VF list --------
GenomeList <- read_excel("../_data/TableS.xlsx", sheet = 2,   skip = 2)
GenomeList <- GenomeList[1:40277,]


VFlist <- read_excel("../_data/TableS.xlsx", sheet = 11, skip = 2)
VFlist[1,3:13]
colnames(VFlist)[3:13] <- as.character(VFlist[1,3:13])
VFlist <- VFlist[-1,]

all( VFlist$`Sequence ID` %in% GenomeList$`Sequence ID a` )


ARGlist <- read_excel("../_data/TableS.xlsx", sheet=5, skip = 2)
as.character(ARGlist[1,4:6])
colnames(ARGlist)[4:6] <- as.character(ARGlist[1,4:6])
ARGlist <- ARGlist[2:453,]

# 整理ARG(- / +)还有MGE(- / +)的表格 ---------

tmp <- 
  merge(
    GenomeList %>% select(`Sequence ID a`),
    ARGlist %>% mutate(ARG = "ARG+") %>% select(`Sequence ID`,ARG),
    by.x = "Sequence ID a", by.y = "Sequence ID",
    all.x = T)

tmp$ARG[is.na(tmp$ARG)] <- "ARG-"


tmp <- 
  merge(
    tmp,
    VFlist %>% mutate(VF="VF+") %>% select(`Sequence ID`, VF),
    by.x = 'Sequence ID a', by.y = "Sequence ID", all.x = T
  )
tmp$VF[is.na(tmp$VF)] <- "VF-"



# 1）chisq.test  ---------
table(tmp$ARG, tmp$VF)
ct <- chisq.test(tmp$ARG, tmp$VF)
ct$statistic
ct$p.value 


# 2） bar plot  ---------
library(ggplot2)

plotDat_phage <- table(tmp$VF, tmp$ARG) %>%
  as.data.frame() %>%
  group_by(Var1) %>%
  mutate(bar.total = sum(Freq)) %>%
  as.data.frame() %>%
  mutate(perc = Freq/bar.total)


ggplot(plotDat_phage %>% filter(Var2 == "ARG+") ) +
  geom_col(aes(x=Var1, y=perc,fill = Var1),width = 0.6, color="black") +
  theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#FFF8BC","#6897C2"))  +
  xlab("") + ylab("ARG-carrying genomes (%)")
ggsave(filename = "Fig6F.VF-ARG_chiSquare_phage.pdf", width = 3, height = 2.2, device = "pdf")
