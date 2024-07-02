
library(data.table)
library(dplyr)
load("../Fig1/Fig1AD.plotSourceDat.RData") # for the Family and colors
color.ncldv.fam <- bind_rows(color.ncldv.fam,
                             cbind.data.frame(family="Overall", color="black"))

library(readxl)
library(data.table)
ARGlist <- read_excel("../_data/TableS.xlsx", sheet = 4, skip= 2)
as.character(ARGlist[1,4:6])
colnames(ARGlist)[4:6] <- as.character(ARGlist[1,4:6])
ARGlist <- ARGlist[2:750,]


# ARG diversity data =======================
tmp <- 
  merge(
    ARGlist %>%
      group_by(`Sequence ID`) %>%
      summarise(n.ARGorf = n()) %>%
      as.data.frame(),
    ARGlist %>%
      select(`Sequence ID`, `Antibiotic type`) %>% unique() %>%
      group_by(`Sequence ID`) %>% 
      summarise(n.ARGtype =n()) %>%
      as.data.frame() ,
    by="Sequence ID"
  )

ARGdiver_df <- merge(ncldv.meta %>% select(`Sequence ID a`, `Sequencing approach`, Family),
                     tmp, by.x = "Sequence ID a", by.y = "Sequence ID", all.x = T)

ARGdiver_df$Family %>% unique()
ARGdiver_df$Family[ARGdiver_df$Family == "n.a."] <- "Unidentified"
ARGdiver_df[is.na(ARGdiver_df)] <- 0



# keep bars with a minimum of 5 ARG-like orfs ==============
Bar.argORFs <- ARGdiver_df %>%
  group_by(`Sequencing approach`, Family) %>%
  summarise(total.argORF = sum(n.ARGorf )) %>%
  mutate(bar =  paste(`Sequencing approach`, Family, sep = "|"))

kept.bars <- Bar.argORFs$bar[Bar.argORFs$total.argORF >= 5] 
kept.bars


# calculate Y limit in y axises  -----
tmp <- ARGdiver_df %>%
  mutate(bar = paste(`Sequencing approach`, Family, sep = "|")) %>%
  filter(bar %in% kept.bars)  %>%
  reshape2::melt(id.vars=c("Sequence ID a", "Sequencing approach","Family", "bar") )%>% 
  group_by(Family, `Sequencing approach`,variable) %>%
  summarise(avg = mean(value), sd = sd(value), 
            n = n()) %>%
  mutate(se = sd/sqrt(n))  
tmp.type <- tmp %>% filter(variable=="n.ARGtype")
maxY.type = max(tmp.type$avg + tmp.type$se)
tmp.orf <- tmp %>% filter(variable=="n.ARGorf")
maxY.orf = max(tmp.orf$avg + tmp.orf$se)




library(ggplot2)
for(tp in c("Isolate","Metagenomic")){
  # tp="Isolate"  
  
  testDat <- ARGdiver_df %>%
    mutate(bar = paste(`Sequencing approach`, Family, sep = "|")) %>%
    filter(bar %in% kept.bars) %>%
    filter(`Sequencing approach` == tp)  
  assign(paste0("testDat.ncldv_",tp), testDat)
  
  
  # statistical grouping of viral families -----
  model_kruskal <- agricolae::kruskal(y=testDat$n.ARGorf, trt= testDat$Family, p.adj = "fdr")
  ARGsubtp.grps <- model_kruskal$groups
  ARGsubtp.grps$n.ARGorf <- sapply(rownames(ARGsubtp.grps), function(x) mean(testDat$n.ARGorf[which(testDat$Family == x)]))
  ARGsubtp.grps <- ARGsubtp.grps %>% arrange(desc(n.ARGorf)) %>% arrange(groups)
  
  model_kruskal <- agricolae::kruskal(y=testDat$n.ARGtype, trt= testDat$Family, p.adj = "fdr")
  ARGtp.grps <- model_kruskal$groups
  ARGtp.grps$n.ARGtype <- sapply(rownames(ARGtp.grps), function(x) mean(testDat$n.ARGtype[which(testDat$Family == x)]))
  ARGtp.grps <- ARGtp.grps %>% arrange(desc(n.ARGtype)) %>% arrange(groups)
  
  
  avg.plotdat <- 
    bind_rows(
      testDat %>%
        group_by(Family) %>%
        summarise(avg.type = mean(n.ARGtype), sd.type = sd(n.ARGtype), 
                  avg.subtype = mean(n.ARGorf), sd.subtype = sd(n.ARGorf),  #是ARGorf，不是subtype
                  n = n()) %>%
        mutate(se.type = sd.type/sqrt(n),
               se.subtype = sd.subtype/sqrt(n)) %>%
        mutate(group.subtype = sapply(Family, function(x)ARGsubtp.grps$groups[which(rownames(ARGsubtp.grps) == x)] ),
               group.type = sapply(Family, function(x)ARGtp.grps$groups[which(rownames(ARGtp.grps) == x)] )) %>%
        as.data.frame(),
      testDat %>%
        group_by(`Sequencing approach`) %>%
        summarise(avg.type = mean(n.ARGtype), sd.type = sd(n.ARGtype), 
                  avg.subtype = mean(n.ARGorf), sd.subtype = sd(n.ARGorf),  #是ARGorf，不是subtype
                  n = n()) %>%
        mutate(se.type = sd.type/sqrt(n),
               se.subtype = sd.subtype/sqrt(n)) %>%
        mutate(Family="Overall",
               group.subtype = NA,
               group.type = NA)
    ) 
  assign(paste0("avg.plotdat.ncldv_",tp), avg.plotdat)
  
  avg.plotdat$Family <- factor(avg.plotdat$Family, levels = c(rownames(ARGtp.grps),"Overall"))
  colors = sapply(c(rownames(ARGtp.grps),"Overall"),
                  function(x) color.ncldv.fam$color[which(color.ncldv.fam$family == x)])
  
  P.tp <- ggplot(avg.plotdat, aes(x=Family, y=avg.type)) +
    geom_errorbar(aes(ymin=avg.type - se.type, ymax=avg.type + se.type), width=0,size=1, color="gray") +
    geom_point(aes(color = Family), size=3 ) +
    theme_bw()+ 
    scale_y_continuous(limits = c(0, maxY.type+0.3)) + 
    scale_color_manual(values = colors) +
    theme_bw() + theme(#axis.text.x =  element_text(angle = 90),
      panel.grid = element_blank()) +
    geom_text(data = avg.plotdat,aes(x=Family, y=maxY.type + 0.2,label=group.type))
  
  P.tp
  assign(paste0("P.type.ncldv_",tp), P.tp)
  
  
}


pdf("Fig2A.No.ARGtype_ncldvs.pdf",  width = 4.8, height = 1.9)
ggarrange(P.type.ncldv_Isolate + theme(legend.position = "none"),
          P.type.ncldv_Metagenomic + theme(legend.position = "none"), widths =  c(0.56, 0.44))
dev.off()


save(kept.bars, ARGdiver_df,  color.ncldv.fam, 
     testDat.ncldv_Isolate, testDat.ncldv_Metagenomic,
     avg.plotdat.ncldv_Isolate, avg.plotdat.ncldv_Metagenomic, 
     P.type.ncldv_Isolate, P.type.ncldv_Metagenomic,
     file = "Fig2A.plots_n.sourceDat_ncldv.RData") 


