

library(dplyr)
load("Fig1AD.plotSourceDat.RData")
color.ncldv.fam <- bind_rows(color.ncldv.fam, 
                             c("family" = "Overall", "color" = "black"))


library(readxl)
library(data.table)
ARGlist <- read_excel("../_data/TableS.xlsx", sheet = 4, skip= 2)
as.character(ARGlist[1,4:6])
colnames(ARGlist)[4:6] <- as.character(ARGlist[1,4:6])
ARGlist <- ARGlist[2:750,]


# calculate percORF.ncldv ----------
percORF.ncldv <- 
  merge(ARGlist %>%
          group_by(`Sequence ID`) %>%
          summarise(numARG=n()) %>%
          as.data.frame(),
        ncldv.meta %>%
          mutate(numORF=`CDS count`,
                 Type=`Sequencing approach`) %>%
          select(`Sequence ID a`,numORF, Type, Family, Habitat),
        by.x="Sequence ID", by.y="Sequence ID a", all=T) %>%
  mutate(perc.orf = numARG/numORF)
percORF.ncldv[is.na(percORF.ncldv)] <- 0
percORF.ncldv$Family[percORF.ncldv$Family == "n.a."] <- "Unidentified"



# plotDat1 for prevalence, plotDat2 for percORF ------------

plotDat1 <- merge(
  ncldv.meta %>%
    group_by(`Sequencing approach`, Family) %>%
    summarise(total.n = n()),
  ncldv.meta %>% filter(`Sequence ID a` %in% ARGlist$`Sequence ID`) %>%
    group_by(`Sequencing approach`, Family) %>%
    summarise(AC.n=n()),
  by=c("Sequencing approach","Family"), all.x = T
) 
plotDat1[is.na(plotDat1)] <- 0
plotDat1$prevalence = plotDat1$AC.n/plotDat1$total.n
plotDat1$Family[plotDat1$Family %in% c("n.a.")] <- "Unidentified"

plotDat1.sub <- 
  bind_rows(
    # Top 5 families with most genomes
    plotDat1 %>%
      group_by(`Sequencing approach`) %>%
      top_n(total.n, n = 5) %>% 
      mutate(bar = paste(`Sequencing approach`, Family, sep = "|")),
    # overall
    plotDat1 %>%
      group_by(`Sequencing approach`) %>%
      summarise(total.n = sum(total.n),
                AC.n = sum(AC.n)) %>%
      mutate(prevalence = AC.n/total.n,
             Family= "Overall") %>%
      mutate(bar = paste(`Sequencing approach`, Family, sep = "|"))
  ) 

kept.bars = plotDat1.sub$bar
unique(plotDat1.sub$Family)


# group types,按照grouping of taxonomy来排序 -----
for(tp in c("Isolate","Metagenomic")){
  # tp = "Metagenomic"
  testDat <- percORF.ncldv %>%
    mutate(bar = paste(Type, Family, sep = "|")) %>%
    filter(bar %in% kept.bars) %>%
    filter(Type == tp)  
  assign(paste0("testDat.ncldv_",tp), testDat)
  
  model_kruskal <- agricolae::kruskal(y=testDat$perc.orf, trt= testDat$Family, p.adj = "fdr")
  groups.dat <- model_kruskal$groups
  
  
  avg.plotdat <- 
    bind_rows( testDat %>%
                 group_by(Family) %>%
                 summarise(avg.per.orf = mean(perc.orf), sd =sd(perc.orf), n=n()) %>%
                 mutate(se=sd/sqrt(n)) %>%
                 mutate(group = sapply(Family,
                                       function(x) groups.dat$groups[which(rownames(groups.dat) == x)])),
               percORF.ncldv %>%
                 group_by(Type) %>%
                 summarise(avg.per.orf = mean(perc.orf), sd=sd(perc.orf),n=n()) %>%
                 mutate(se=sd/sqrt(n), 
                        Family="Overall") %>%
                 filter(Type == tp))
  
  lvls.df <- avg.plotdat %>% arrange(desc(avg.per.orf)) %>% arrange(group)
  
  habt.lvl <- lvls.df$Family
  avg.plotdat$Family <- factor(avg.plotdat$Family, levels = habt.lvl)
  assign(paste0("plotDat2.ncldv_",tp), avg.plotdat)
  
  colors = sapply(habt.lvl, function(x) color.ncldv.fam$color[which(color.ncldv.fam$family == x)])
  
  
  P <- ggplot(avg.plotdat) +
    geom_col(aes(x=Family, y=avg.per.orf, fill=Family), width = 0.8 ) +
    geom_errorbar(aes(x=Family, y=avg.per.orf, ymin=avg.per.orf-se, ymax=avg.per.orf+se), width=.1, color="gray") +
    theme_bw()+ 
    #scale_y_continuous(labels = scales::percent,limits = c(0,0.007)) + # 先不设limits，看了两个type的图后统一: isolates: c(0,0.007)
    scale_fill_manual(values = colors) +
    theme_bw() + theme(panel.grid = element_blank()) +
    geom_text(aes(x=Family, y=0.0065, label=group)) # 根据前面的ylimits设置y值
  
  P
  assign(paste0("P2.ncldv_",tp), P)
  
  
  # 按新的排序把p1的图也重新画一下
  plotDat1.ncldv <- plotDat1.sub %>% filter( `Sequencing approach` == tp)
  plotDat1.ncldv$Family <- factor(plotDat1.ncldv$Family, levels = habt.lvl)
  assign(paste0("plotDat1.ncldv_",tp), plotDat1.ncldv)
  
  P1 <- ggplot(plotDat1.ncldv) +
    geom_col(aes(x=Family, y=prevalence , fill=Family),  width = 0.8) +
    theme_bw()+ 
    scale_y_continuous(labels = scales::percent,limits = c(0,1)) + # 先不设limits，看了两个type的图后统一设: isolates: 0-1
    scale_fill_manual(values = colors) +
    theme_bw() + theme(panel.grid = element_blank()) 
  P1
  
  assign(paste0("P1.nlcdv_",tp), P1)
  
  
}


# rearrange x axis (the viral families ) according to genome number ===================

plotDat1.ncldv_Isolate <- plotDat1.ncldv_Isolate %>% arrange(desc(total.n))
IG.lvls <- c(as.character(plotDat1.ncldv_Isolate$Family[plotDat1.ncldv_Isolate$Family != "Overall"]), "Overall")
plotDat1.ncldv_Isolate$Family <- factor(plotDat1.ncldv_Isolate$Family, 
                                        levels = IG.lvls)
colors <- sapply(IG.lvls, function(x) color.ncldv.fam$color[which(color.ncldv.fam$family == x)])
P1_IG <- ggplot(plotDat1.ncldv_Isolate) +
  geom_col(aes(x=Family, y=prevalence , fill=Family),  width = 0.8) +
  theme_bw()+ 
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) + # 先不设limits，看了两个type的图后统一设: isolates: 0-1
  scale_fill_manual(values = colors) +
  theme_bw() + theme(panel.grid = element_blank()) 
P1_IG


plotDat2.ncldv_Isolate$Family <- factor(plotDat2.ncldv_Isolate$Family, 
                                        levels = IG.lvls)
P2_IG <- ggplot(plotDat2.ncldv_Isolate) +
  geom_col(aes(x=Family, y=avg.per.orf, fill=Family), width = 0.8 ) +
  geom_errorbar(aes(x=Family, y=avg.per.orf, ymin=avg.per.orf-se, ymax=avg.per.orf+se), width=.1, color="gray") +
  theme_bw()+ 
  scale_y_continuous(labels = scales::percent,limits = c(0,0.007)) + # 先不设limits，看了两个type的图后统一: isolates: c(0,0.007)
  scale_fill_manual(values = colors) +
  theme_bw() + theme(panel.grid = element_blank()) +
  geom_text(aes(x=Family, y=0.0065, label=group)) # 根据前面的ylimits设置y值

P2_IG




plotDat1.ncldv_Metagenomic <- plotDat1.ncldv_Metagenomic %>% arrange(desc(total.n))
MAG.lvls <- c(as.character(plotDat1.ncldv_Metagenomic$Family[plotDat1.ncldv_Metagenomic$Family != "Overall"]), "Overall")
plotDat1.ncldv_Metagenomic$Family <- factor(plotDat1.ncldv_Metagenomic$Family, 
                                            levels = MAG.lvls)
colors <- sapply(MAG.lvls, function(x) color.ncldv.fam$color[which(color.ncldv.fam$family == x)])
P1_MAG <- ggplot(plotDat1.ncldv_Metagenomic) +
  geom_col(aes(x=Family, y=prevalence , fill=Family),  width = 0.8) +
  theme_bw()+ 
  scale_y_continuous(labels = scales::percent,limits = c(0,1)) + # 先不设limits，看了两个type的图后统一设: isolates: 0-1
  scale_fill_manual(values = colors) +
  theme_bw() + theme(panel.grid = element_blank()) 
P1_MAG

plotDat2.ncldv_Metagenomic$Family <- factor(plotDat2.ncldv_Metagenomic$Family,
                                            levels = MAG.lvls)
P2_MAG <- ggplot(plotDat2.ncldv_Metagenomic) +
  geom_col(aes(x=Family, y=avg.per.orf, fill=Family), width = 0.8 ) +
  geom_errorbar(aes(x=Family, y=avg.per.orf, ymin=avg.per.orf-se, ymax=avg.per.orf+se), width=.1, color="gray") +
  theme_bw()+ 
  scale_y_continuous(labels = scales::percent,limits = c(0,0.007)) + # 先不设limits，看了两个type的图后统一: isolates: c(0,0.007)
  scale_fill_manual(values = colors) +
  theme_bw() + theme(panel.grid = element_blank()) +
  geom_text(aes(x=Family, y=0.0065, label=group)) # 根据前面的ylimits设置y值

P2_MAG


pdf("Fig1BC.ncldv_top5_ARGcontent.pdf",  width = 8, height = 3)
ggarrange(ggarrange(P1_IG + theme(legend.position = "none"),
                    P1_MAG + theme(legend.position = "none"), widths = c(0.5, 0.5)),
          ggarrange(P2_IG + theme(legend.position = "none"),
                    P2_MAG+ theme(legend.position = "none"), widths = c(0.5, 0.5)),
          widths = c(0.5,0.5))
dev.off()


