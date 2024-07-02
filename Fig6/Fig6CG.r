
library(readxl)
library(data.table)
library(dplyr)
library(scales)
library(ggplot2)


# NCLDV =============================================
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


plotDat.ncldv <- GenomeList %>%
  filter(`Sequence ID a` %in% ARGlist$`Sequence ID`) %>%
  filter(`Sequence ID a` %in% VFlist$`Sequence ID`) %>%
  group_by(Family) %>%
  summarise(count=n()) %>%
  as.data.frame() %>%
  arrange(Family) %>%
  mutate(fraction = count/sum(count))  %>%
  mutate(ymax= cumsum(fraction)) %>%
  mutate(ymin= c(0, head(ymax, n=-1))) %>%
  mutate(labelPosition = (ymax + ymin) / 2) %>%
  mutate(perc = percent(fraction)) %>%
  mutate(label = paste0(Family, " (", perc, ")"))
plotDat.ncldv$Family[plotDat.ncldv$Family == "n.a."] <- "Unidentified"


load("../Fig1/Fig1AD.plotSourceDat.RData") # for the colors
colors.ncldv <- setNames(color.ncldv.fam$color, color.ncldv.fam$family)
ggplot(plotDat.ncldv, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=1, fill=Family)) + # xmin ~ xmax is where the ring is 
  geom_rect(color="white",lwd=0.5) +
  # geom_text( x=3.5, aes(y=labelPosition, label=perc), color="white", size=3) + # x here controls label position (inner / outer)
  geom_text( x=5, aes(y=labelPosition, label=label), color="black", size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = colors.ncldv) +
  coord_polar(theta="y") +
  xlim(c(1, 5)) + # 可调节donut thickness
  theme_void() +
  theme(legend.position = "none")
ggsave(filename = "Fig6C.ARG-VF_ncldvFamilies.pdf", device = "pdf", width = 3, height = 3)





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

plotDat.phage <- GenomeList %>%
  filter(`Sequence ID a` %in% ARGlist$`Sequence ID`) %>%
  filter(`Sequence ID a` %in% VFlist$`Sequence ID`) %>%
  group_by(Family) %>%
  summarise(count=n()) %>%
  as.data.frame() %>%
  arrange(Family) %>%
  mutate(fraction = count/sum(count))  %>%
  mutate(ymax= cumsum(fraction)) %>%
  mutate(ymin= c(0, head(ymax, n=-1))) %>%
  mutate(labelPosition = (ymax + ymin) / 2) %>%
  mutate(perc = percent(fraction)) %>%
  mutate(label = paste0(Family, " (", perc, ")"))
plotDat.phage$Family[plotDat.phage$Family == "n.a."] <- "Unidentified"


colors.phage <- setNames(c( "darkorange3", "#EAB33AFF","#B1746FFF","lightgray"  ), 
                         c("Sphaerolipoviridae","Herelleviridae","Straboviridae","Unidentified"))

ggplot(plotDat.phage, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=1, fill=Family)) + # xmin ~ xmax is where the ring is 
  geom_rect(color="white",lwd=0.5) +
  # geom_text( x=3.5, aes(y=labelPosition, label=perc), color="white", size=3) + # x here controls label position (inner / outer)
  geom_text( x=5, aes(y=labelPosition, label=label), color="black", size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = colors.phage) +
  coord_polar(theta="y") +
  xlim(c(1, 5)) + # 可调节donut thickness
  theme_void() +
  theme(legend.position = "none")
ggsave(filename = "Fig6G.ARG-VF_phageFamilies.pdf", device = "pdf", width = 3, height = 3)
