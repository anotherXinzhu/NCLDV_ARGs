

library(data.table)
library(readxl)
library(dplyr)
library(scales)



# Viral genome list ===========================================
ncldv.meta <- read_excel("../_data/TableS.xlsx", sheet = 1, skip = 2)
ncldv.meta <- ncldv.meta[1:1416,] 

phage.meta <- read_excel("../_data/TableS.xlsx", sheet = 2, skip = 2) 
phage.meta <- phage.meta[1:40277,] 


# assign colors to each taxonomy ===========================================
library(paletteer) 
unique(ncldv.meta$Family)
color.ncldv.fam <- 
  cbind.data.frame(
    family= c("Asfarviridae","Iridoviridae","Marseilleviridae","Mimiviridae",
               "Pandoraviridae","Phycodnaviridae","Pithoviridae","Poxviridae",
               "Prasinoviridae","Unidentified","Others",
               "Mininucleoviridae", "Coccolithoviridae"),
    color = c("#042F49FF", "#2EA7E0", "#3394BBFF", "#8DB576FF",
              "#418877FF", "#146252FF", "#B06429FF", "#E5AE8DFF",
              "#F6EBD1FF","lightgray","darkgray", 
              "#920783","#C0392BFF"),
    stringsAsFactors=F
  ) 

# top 10 phage families with most genomes-----
phage2show <- names(table(phage.meta$Family)[order(table(phage.meta$Family), decreasing = T)][1:11])   

table(phage.meta$Family)[phage2show]
color.phage.fam <- 
  cbind.data.frame(
    # family= c("Herelleviridae","Herpesviridae","Myoviridae","Unidentified","Siphoviridae","Microviridae","Inoviridae" ,     
    #            "Podoviridae","Papillomaviridae", "Circoviridae","Flaviviridae","Reoviridae","Retroviridae","Others","n.a."),
    family= c("Autographiviridae","Genomoviridae","Anelloviridae","Unidentified","Sedoreoviridae","Microviridae","Inoviridae" ,     
               "Podoviridae","Papillomaviridae", "Circoviridae","Flaviviridae","Schitoviridae","Retroviridae","Others","n.a."),
    color = c("darkorange3", "#EAB33AFF", "#B1746FFF", "lightgray", "#8A8B79FF", "#EB9FAFFF", "#FFEFA9FF",
              "#AD2E37FF","#882B6AFF","#E61673","gold","#6A5188FF","#4882B6FF","darkgray","darkgray"),
    stringsAsFactors=F
  ) # top 10 num sequence or ARG-carrying genome >= 5

phage2show[!phage2show %in% color.phage.fam$family]
color.phage.fam$family[!color.phage.fam$family %in% phage2show]

color.phage.fam$color[which(color.phage.fam$family %in% phage2show)] %>% show_col()




# plotting =======================================================
# ncldv -----------------------

tmpPD.ncldv <- ncldv.meta 
tmpPD.ncldv$Family %>% unique
tmpPD.ncldv$Family[tmpPD.ncldv$Family %in% c("-","n.a.")] <- "Unidentified"

# MAG
plotDat.ncldv_MAG <-
  tmpPD.ncldv %>%
  filter(`Sequencing approach` == "Metagenomic") %>%
  group_by(Family) %>% 
  summarise(count = n()) %>%
  as.data.frame() %>%
  arrange(Family) %>%
  mutate(fraction = count/sum(count)) %>%
  mutate(ymax= cumsum(fraction)) %>%
  mutate(ymin= c(0, head(ymax, n=-1))) %>%
  mutate(labelPosition = (ymax + ymin) / 2) %>%
  mutate(perc = percent(fraction)) %>%
  mutate(label = paste0(Family, " (", perc, ")"))


library(ggplot2)
colors <- sapply(plotDat.ncldv_MAG$Family,
                 function(x) color.ncldv.fam$color[which(color.ncldv.fam$family == x)])
p_ncldv.MAG <- ggplot(plotDat.ncldv_MAG, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Family)) + # xmin ~ xmax is where the ring is 
  geom_rect(color="white",lwd=0.5) +
  geom_text( x=5, aes(y=labelPosition, label=label), color="black", size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = colors) +
  coord_polar(theta="y") +
  xlim(c(1, 5)) + # adjust donut thickness
  theme_void() +
  theme(legend.position = "none")

p_ncldv.MAG


# IGs
plotDat.ncldv_IG <-
  tmpPD.ncldv %>%
  filter(`Sequencing approach` == "Isolate") %>%
  group_by(Family) %>% 
  summarise(count = n()) %>%
  as.data.frame() %>%
  arrange(Family) %>%
  mutate(fraction = count/sum(count)) %>%
  mutate(ymax= cumsum(fraction)) %>%
  mutate(ymin= c(0, head(ymax, n=-1))) %>%
  mutate(labelPosition = (ymax + ymin) / 2) %>%
  mutate(perc = percent(fraction)) %>%
  mutate(label = paste0(Family, " (", perc, ")"))

colors <- sapply(plotDat.ncldv_IG$Family,
                 function(x) color.ncldv.fam$color[which(color.ncldv.fam$family == x)])
p_ncldv.IG <- ggplot(plotDat.ncldv_IG, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Family)) + # xmin ~ xmax is where the ring is 
  geom_rect(color="white",lwd=0.5) +
  geom_text( x=5, aes(y=labelPosition, label=label), color="black", size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = colors) +
  coord_polar(theta="y") +
  xlim(c(1, 5)) + # adjust donut thickness
  theme_void() +
  theme(legend.position = "none")


library(ggpubr)
ggarrange(p_ncldv.IG, p_ncldv.MAG)





# phage -----------------------
tmpPD.phage <- phage.meta 
unique(tmpPD.phage$Family)
tmpPD.phage$Family[tmpPD.phage$Family %in% c("-","n.a.")] <- "Unidentified"
tmpPD.phage$Family[!tmpPD.phage$Family %in% c(phage2show, "Unidentified")] <- "Others"


# MAG
plotDat.phage_MAG <-
  tmpPD.phage %>%
  filter(`Sequencing approach` == "Metagenomic") %>%
  group_by(Family) %>% 
  summarise(count = n()) %>%
  as.data.frame() %>%
  arrange(Family) %>%
  mutate(fraction = count/sum(count)) %>%
  mutate(ymax= cumsum(fraction)) %>%
  mutate(ymin= c(0, head(ymax, n=-1))) %>%
  mutate(labelPosition = (ymax + ymin) / 2) %>%
  mutate(perc = percent(fraction)) %>%
  mutate(label = paste0(Family, " (", perc, ")"))

library(ggplot2)
colors <- sapply(plotDat.phage_MAG$Family,
                 function(x) color.phage.fam$color[which(color.phage.fam$family == x)])
p_phage.MAG <- ggplot(plotDat.phage_MAG, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Family)) + # xmin ~ xmax is where the ring is 
  geom_rect(color="white",lwd=0.5) +
  geom_text( x=5, aes(y=labelPosition, label=label), color="black", size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = colors) +
  coord_polar(theta="y") +
  xlim(c(1, 5)) + # adjust donut thickness
  theme_void() +
  theme(legend.position = "none")

p_phage.MAG



# IG
plotDat.phage_IG <-
  tmpPD.phage %>%
  filter(`Sequencing approach` == "Isolate") %>%
  group_by(Family) %>% 
  summarise(count = n()) %>%
  as.data.frame() %>%
  arrange(Family) %>%
  mutate(fraction = count/sum(count)) %>%
  mutate(ymax= cumsum(fraction)) %>%
  mutate(ymin= c(0, head(ymax, n=-1))) %>%
  mutate(labelPosition = (ymax + ymin) / 2) %>%
  mutate(perc = percent(fraction)) %>%
  mutate(label = paste0(Family, " (", perc, ")"))


library(ggplot2)
colors <- sapply(plotDat.phage_IG$Family,
                 function(x) color.phage.fam$color[which(color.phage.fam$family == x)])
show_col(colors)
p_phage.IG <- ggplot(plotDat.phage_IG, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Family)) + # xmin ~ xmax is where the ring is 
  geom_rect(color="white",lwd=0.5) +
  geom_text( x=5, aes(y=labelPosition, label=label), color="black", size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = colors) +
  coord_polar(theta="y") +
  xlim(c(1, 5)) + # adjust donut thickness
  theme_void() +
  theme(legend.position = "none")

p_phage.IG



library(ggpubr)
ggarrange(p_phage.IG, p_phage.MAG)


# integrate plots  -----
pdf("Fig1AD.num.genomes.pdf", width = 7, height = 5.5)
ggarrange(
  ggarrange(p_ncldv.IG + ggtitle("NCLDV", subtitle = "IG"), p_ncldv.MAG + ggtitle("NCLDV", subtitle = "MAG")),
  ggarrange(p_phage.IG + ggtitle("Phage", subtitle = "IG"), p_phage.MAG + ggtitle("Phage", subtitle = "MAG")),
  nrow = 2
)
dev.off()




# legend NCLDv -----
colors <- sapply(color.ncldv.fam$family,
                 function(x) color.ncldv.fam$color[which(color.ncldv.fam$family == x)])

ncldv.leg <- 
  ggplot(color.ncldv.fam) +
  geom_tile(aes(x=family, y=0, fill=family), shape=21, size=3) +
  theme_void() +
  theme(legend.position="bottom") + guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = colors)
ncldv.leg

pdf("Fig1A.ncldv-legend.pdf", width = 14, height = 2)
ncldv.leg
dev.off()


# legend phage -----
colors <- sapply(c(phage2show[phage2show != "n.a."],"Unidentified"),
                 function(x) color.phage.fam$color[which(color.phage.fam$family == x)])

phage.leg <- 
  ggplot(color.phage.fam) +
  geom_tile(aes(x=family, y=0, fill=family)) +
  theme_void() +
  theme(legend.position="bottom") + guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = colors)
phage.leg

pdf("Fig1D.phage-legend.pdf", width = 14, height = 2)
phage.leg
dev.off()



save(phage.meta, ncldv.meta,  
     color.ncldv.fam,color.phage.fam, 
     file = "Fig1AD.plotSourceDat.RData")


