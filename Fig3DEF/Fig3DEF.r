library(data.table)


# Fig3D ===============================================
dat_dfr <- fread("dfr_dNdS.txt",data.table = F)

library(ggplot2)
library(gghalves)
library(ggpubr)

myList <- combn(unique(dat_dfr$kingdom), 2) %>% as.data.frame(stringsAsFactors=F) %>% as.list()
?stat_compare_means
ggboxplot(dat_dfr, x="kingdom",y="Ka/Ks") +
  stat_compare_means(comparisons = myList) +
  theme_bw() # check significance


# manually assign statistical groups 
dat_dfr$kingdom <- factor(dat_dfr$kingdom, levels = c("ncldv","phage","bacteria","eukaryote"))

# assign groups, then calculate wilcox
dat_dfr$group <- sapply(dat_dfr$kingdom,
                         function(x) if(x == "ncldv") "c" else if(x %in% c("bacteria","phage")) "b" else "a")
myList <- combn(unique(dat_dfr$group), 2) %>% as.data.frame(stringsAsFactors=F) %>% as.list()
for(gs in myList){
  # print(gs)
  
  testD <- dat_dfr %>% filter(group %in% gs)
  w = wilcox.test(testD$`Ka/Ks`~testD$group)
  print(paste0(paste(gs, collapse = " vs. "),  ", wilcox.p: ", w$p.value))
}

groups.dat <- dat_dfr %>% select(kingdom, group) %>% unique
colors <- c("ncldv" = "#D8695B", "phage"="#C1B3CC",
            "eukaryote" = "#ECD494", "bacteria" = "#58855C")

ggplot(dat_dfr, aes(kingdom, `Ka/Ks`, fill=kingdom)) + 
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                   side='R',adjust=1.2,trim=F,color="white",alpha=0.8) +
  theme_bw() +
  theme(legend.position = "none") + 
  geom_text(data=groups.dat, aes(x=kingdom, y=0.6, label=group)) +
  geom_boxplot( width = .15,  outlier.shape = NA )    +
  scale_fill_manual(values = colors) +
  ylab("Ka/Ks of dfr")
ggsave("Fig3D.dNdS_dfr.pdf", width = 3.3, height = 3.2, device = "pdf")


# Fig3E ===============================================

dat_ABCF <- fread("ABCF_dNdS.txt",data.table = F)

myList <- combn(unique(dat_ABCF$kingdom), 2) %>% as.data.frame(stringsAsFactors=F) %>% as.list()
ggboxplot(dat_ABCF %>% filter(`Ka/Ks` < 5),
          x="kingdom",y="Ka/Ks") +
  stat_compare_means(comparisons = myList) +
  theme_bw()

# assign groups then calculate wilxocon p
dat_ABCF$kingdom <- factor(dat_ABCF$kingdom, levels = c("eukaryotes","ncldv","phage","bacteria"))
dat_ABCF$group <- sapply(dat_ABCF$kingdom,
                         function(x){
                           if(x == "bacteria") {
                             "a" 
                           } else if(x %in% c("phage")) {
                             "b"
                           } else if(x %in% "ncldv"){
                             "c"
                           }else "d"
                         } )
myList <- combn(unique(dat_ABCF$group), 2) %>% as.data.frame(stringsAsFactors=F) %>% as.list()
for(gs in myList){
  # print(gs)
  
  testD <- dat_ABCF %>% filter(group %in% gs)
  w = wilcox.test(testD$`Ka/Ks`~testD$group)
  print(paste0(paste(gs, collapse = " vs. "),  ", wilcox.p: ", w$p.value))
}


groups.dat_ABCF <- dat_ABCF %>% select(kingdom, group) %>% unique
colors <- c("ncldv" = "#D8695B", "phage"="#C1B3CC",
            "eukaryotes" = "#ECD494", "bacteria" = "#58855C")

p<- ggplot(dat_ABCF %>% filter(`Ka/Ks` < 1.5), 
           aes(kingdom, `Ka/Ks`, fill=kingdom)) + 
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                   side='R',adjust=1.2,trim=F,color="white",alpha=0.8,scale = "width") +
  # geom_flat_violin(position = position_nudge(x = .2, y = 0)) + 
  # geom_jitter(alpha = 0.5,width = 0.15) + 
  theme_bw() +
  theme(legend.position = "none") + 
  geom_text(data=groups.dat_ABCF, aes(x=kingdom, y=1.3, label=group)) +
  geom_boxplot( width = .1,  outlier.shape = NA )   +
  coord_cartesian(ylim = c(-0.1, 1.5)) +
  scale_fill_manual(values = colors) +
  ylab("Ka/Ks of ABC-F")
p
ggsave("Fig3E.dNdS_ABCF.pdf", width = 3.3, height = 3.2, device = "pdf")



# Fig 3F ==================================

dat_ileS <- fread("ileS_dNdS.txt", data.table = F)

# check statistical significance 
myList <- combn(unique(dat_ileS$kingdom), 2) %>% as.data.frame(stringsAsFactors=F) %>% as.list()
?stat_compare_means
ggboxplot(dat_ileS %>% filter(`Ka/Ks` < 5),
          x="kingdom",y="Ka/Ks") +
  stat_compare_means(comparisons = myList) +
  theme_bw()


# assign groups then calculate wilxocon p
dat_ileS$kingdom <- factor(dat_ileS$kingdom, levels = c("ncldv","eukaryotes","bacteria"))
dat_ileS$group <- sapply(dat_ileS$kingdom,
                         function(x) if(x == "ncldv") "c" else if(x %in% c("eukaryotes")) "b" else "a")
myList <- combn(unique(dat_ileS$group), 2) %>% as.data.frame(stringsAsFactors=F) %>% as.list()
for(gs in myList){
  # print(gs)
  
  testD <- dat_ileS %>% filter(group %in% gs)
  w = wilcox.test(testD$`Ka/Ks`~testD$group)
  print(paste0(paste(gs, collapse = " vs. "),  ", wilcox.p: ", w$p.value))
}

groups.dat_ileS <- dat_ileS %>% select(kingdom, group) %>% unique
colors <- c("ncldv" = "#D8695B", "eukaryotes" = "#ECD494", "bacteria" = "#58855C")

p <- ggplot(dat_ileS ,
            aes(kingdom, `Ka/Ks`, fill=kingdom)) + 
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                   side='R',adjust=1.2,trim=F,color="white",alpha=0.8, scale = "width") +
  theme_bw() +
  theme(legend.position = "none") + 
  geom_text(data=groups.dat_ileS, aes(x=kingdom, y=0.9, label=group)) +
  geom_boxplot( width = .10,  outlier.shape = NA )    +
  scale_fill_manual(values = colors) +
  ylab("Ka/Ks of ileS")
p
ggsave("Fig3E.dNdS_ileS.pdf", width = 3.3, height = 3.2, device = "pdf")
