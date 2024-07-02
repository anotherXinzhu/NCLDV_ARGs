
library(readxl)
library(dplyr)
library(ggplot2)


# genome list ===========
GenomeList <- read_excel("../_data/TableS.xlsx", sheet = 2, skip = 2)
GenomeList <- GenomeList[1:40277,]


# Fig 5E. MGE+ genomes ================================
MGE.list <- fread("MGElist_phage.txt", data.table = F)
tmp <- merge(
  GenomeList %>% select(`Sequence ID a`),
  MGE.list %>%
    mutate(available = "MGE+") %>%
    select(genome, available),
  by.x="Sequence ID a", 
  by.y="genome", 
  all=T
) %>% unique

tmp$available[is.na(tmp$available) ] <- "MGE-"

plotDat.5E <- 
  tmp %>%
  group_by(available) %>%
  summarise(n = n()) %>%
  mutate(sum = sum(n)) %>%
  mutate(perc = n/sum)

ggplot(plotDat.5E) +
  geom_col(aes(x=available, y=perc, fill=available), width = 0.6, color="black") +
  scale_fill_manual(values = c("#FFF8BC","#6897C2")) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  ylab("Genomes (%)") + xlab("")
ggsave("Fig5E.MGEPrevalence_phage.pdf", device = "pdf", width = 3, height = 2.2)


# Fig 5F ARG MGE Chi-squared =====================================
ARGlist <- read_excel("../_data/TableS.xlsx", sheet=5, skip = 2)
as.character(ARGlist[1,4:6])
colnames(ARGlist)[4:6] <- as.character(ARGlist[1,4:6])
ARGlist <- ARGlist[2:454,]

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
    MGE.list %>% mutate(MGE="MGE+") %>% select(genome, MGE),
    by.x = 'Sequence ID a', by.y = "genome", all.x = T
  )
tmp$MGE[is.na(tmp$MGE)] <- "MGE-"

# 1）chisq.test  ---------
table(tmp$ARG, tmp$MGE)
ct <- chisq.test(tmp$ARG, tmp$MGE)# 
ct$statistic
ct$p.value  

# 2） bar plot  ---------
plotDat.5F <- table(tmp$MGE, tmp$ARG) %>%
  as.data.frame() %>%
  group_by(Var1) %>%
  mutate(bar.total = sum(Freq)) %>%
  as.data.frame() %>%
  mutate(perc = Freq/bar.total)

ggplot(plotDat.5F %>% filter(Var2 == "ARG+"),  
             aes(x=Var1, y=perc)) +
  geom_col(aes(fill = Var1),width = 0.6, color="black") +
  theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#FFF8BC","#6897C2"))  +
  xlab("") + ylab("ARG-carrying genomes (%)")
ggsave("Fig5F.MGE-ARG_chiSquare_phage.pdf", device = "pdf", width = 3, height = 2.2)



# Fig 5G. difference of ARG% per CDS between MGE- and MGE+ genomes =========================
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
  mutate(MGE = sapply(`Sequence ID a`, function(x) if(x %in% MGE.list$genome) "MGE+" else "MGE-"))

w = wilcox.test(ARGNo_df$ARGperc ~ ARGNo_df$MGE)

plotDat.5G <- ARGNo_df %>%
  group_by(MGE) %>%
  summarise(avg=mean(ARGperc), sd=sd(ARGperc), n=n()) %>%
  as.data.frame() %>%
  mutate(se=sd/sqrt(n))

ggplot(plotDat.5G, aes(x=MGE, y=avg)) +
  geom_col(aes(fill=MGE), width = 0.6, color="black") +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.1, color="gray") +
  scale_y_continuous(labels = scales::percent,
                     #limits = c(0,0.007),   # 看了两个type的图后统一: isolates: c(0,0.007)
  ) +
  scale_fill_manual(values = c("#FFF8BC","#6897C2")) +
  theme_classic() +
  xlab("") + ylab("ARG-like ORFs per genome (%)") 
ggsave("Fig5G.MGE-ARGperc_wilcox_phage.pdf", device = "pdf", width = 3, height = 2.2)


# Fig 5H. ARG-MGE distances ==========================================
ARG.MGE_Cooccurrence <- fread("ARG.MGE_Cooccurrence_phage.txt", data.table = F)

plotDat.5H <- ARG.MGE_Cooccurrence %>%
  dplyr::group_by(ARGorf) %>%
  dplyr::summarise(MGE = Element[which(`ARG-MGE.distance` ==  min(`ARG-MGE.distance`))],
                   MGE.type = MG.type[which(`ARG-MGE.distance` ==  min(`ARG-MGE.distance`))],
                   closestMGE = min(`ARG-MGE.distance`),
                   ARG.MGE.relation = `ARG-MGE.relation`[which(`ARG-MGE.distance` ==  min(`ARG-MGE.distance`))]) %>%
  as.data.frame() %>% unique

ggplot(plotDat.5H, aes(closestMGE)) +
  geom_histogram(binwidth = 5000, fill = "skyblue", color = "gray", size=0.2) +
  scale_x_continuous(labels = function(x) x/1000, # x轴按照1000bp的单位显示
                     breaks = seq(0, 200000, by = 50000) # 设置x轴的标签间隔
  ) +
  coord_cartesian(xlim = c(0, 200000)) +  # set the axis range without changing values
  labs(x = "ARG-MGE distance（kb）", y="Count of ARG-MGE pairs") + 
  theme_classic() #+ theme(panel.grid = element_blank())
ggsave("Fig5H.ARG-MGEdist_phage.pdf", device = 'pdf', width = 3, height = 2.2)

