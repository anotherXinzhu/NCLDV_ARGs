
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


# calculate ARG perc --------------------

ARG.no <- ARGlist %>% 
  group_by(`Sequence ID`) %>%
  summarise(No.ARG = n()) %>%
  as.data.frame()

ARGNo_df <- merge(GenomeList %>% select(`Sequence ID a`, `CDS count`), 
                  ARG.no, 
                  by.x = "Sequence ID a", by.y = "Sequence ID",
                  all.x = T)
ARGNo_df$No.ARG[is.na(ARGNo_df$No.ARG)] <- 0

ARGNo_df <- ARGNo_df %>%
  mutate(ARGperc=No.ARG/`CDS count`) %>% 
  mutate(VF = sapply(`Sequence ID a`,
                     function(x) if(x %in% VFlist$`Sequence ID`) "VF+" else "VF-"))

w = wilcox.test(ARGNo_df$ARGperc ~ ARGNo_df$VF)

plotDat_ncldv <- ARGNo_df %>%
  group_by(VF) %>%
  summarise(avg=mean(ARGperc), sd=sd(ARGperc), n=n()) %>%
  as.data.frame() %>%
  mutate(se=sd/sqrt(n))

ggplot(plotDat_ncldv, aes(x=VF, y=avg)) +
  geom_col(aes(fill=VF), width = 0.6, color="black") +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.1, color="gray") +
  scale_y_continuous(labels = scales::percent,
                     #limits = c(0,0.007),   # 看了两个type的图后统一: isolates: c(0,0.007)
  ) +
  scale_fill_manual(values = c("#FFF8BC","#6897C2")) +
  theme_classic() +
  xlab("") + ylab("ARG-like ORFs per genome (%)") 
ggsave("Fig6D.VF-ARGperc_wilcox_ncldv.pdf",device = "pdf", width = 3, height = 2.2)


# phages =======================================

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



ARG.no <- ARGlist %>% 
  group_by(`Sequence ID`) %>%
  summarise(No.ARG = n()) %>%
  as.data.frame()

ARGNo_df <- merge(GenomeList %>% select(`Sequence ID a`, `CDS count`), 
                  ARG.no, 
                  by.x = "Sequence ID a", by.y = "Sequence ID",
                  all.x = T)
ARGNo_df$No.ARG[is.na(ARGNo_df$No.ARG)] <- 0
ARGNo_df <- ARGNo_df %>% 
  mutate(ARGperc=No.ARG/`CDS count`) %>% 
  mutate(VF = sapply(`Sequence ID a`,
                     function(x) if(x %in% VFlist$`Sequence ID`) "VF+" else "VF-"))

w = wilcox.test(ARGNo_df$ARGperc ~ ARGNo_df$VF)
w

plotDat_phage <- ARGNo_df %>%
  group_by(VF) %>%
  summarise(avg=mean(ARGperc), sd=sd(ARGperc), n=n()) %>%
  as.data.frame() %>%
  mutate(se=sd/sqrt(n))

ggplot(plotDat_phage, aes(x=VF, y=avg)) +
  geom_col(aes(fill=VF), width = 0.6, color="black") +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.1, color="gray") +
  scale_y_continuous(labels = scales::percent,
                     #limits = c(0,0.007),   # 看了两个type的图后统一: isolates: c(0,0.007)
  ) +
  scale_fill_manual(values = c("#FFF8BC","#6897C2")) +
  theme_classic() +
  xlab("") + ylab("ARG-like ORFs per genome (%)") 
ggsave("Fig6H.VF-ARGperc_wilcox_phage.pdf",device = "pdf", width = 3, height = 2.2)
