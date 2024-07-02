
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(scales)
library(RColorBrewer )

# Colors ------ 
library(paletteer) 
paletteer_d("tidyquant::tq_light")
colourCount = 15
getPalette = colorRampPalette(paletteer_d("tidyquant::tq_light"))
show_col(getPalette(colourCount)) 

colors.ARGtype <- 
  cbind.data.frame(
    ARGtype = c( "aminoglycoside","aminocoumarin","phenicol","multidrug","polymyxin","tetracycline","trimethoprim","mupirocin", 
                 "sulfonamide","peptide","MLS","rifamycin","disinfecting agents and antiseptics","quinolone","nitroimidazole","pleuromutilin","Others"),
    color = c("#2C3E90",getPalette(colourCount), "gray"),
    stringsAsFactors=F
  )


# the family levels (to present in figures) -----
load("../Fig1/Fig1AD.plotSourceDat.RData") # for the meta data
load("Fig2B.plots_n.sourceDat_phage.RData") # this is for the levels of kept bars and family levels
fam.lvl_Isolate <- (avg.plotdat.phage_Isolate %>% arrange(desc(avg.type)) %>% arrange(group.type))$Family
fam.lvl_Metagenomic <- (avg.plotdat.phage_Metagenomic %>% arrange(desc(avg.type)) %>% arrange(group.type))$Family
remove(ARGdiver_df)
remove(list = ls()[grepl("^P", ls())])
remove(list = ls()[grepl("^testDat", ls())])
remove(list = ls()[grepl("^avg.plotdat.ncldv", ls())])


# the ARG list ----
library(readxl)

ARGlist <- read_excel("../_data/TableS.xlsx",sheet = 5, skip=2)
as.character(ARGlist[1,4:6])
colnames(ARGlist)[4:6] <- as.character(ARGlist[1,4:6])
ARGlist <- ARGlist[2:454,]

ARGlist <- ARGlist %>%
  mutate(Family = sapply(`Sequence ID`, function(x) phage.meta$Family[which(phage.meta$`Sequence ID a`== x)])) %>%
  mutate(`Sequencing approach` = sapply( `Sequence ID`, function(x)phage.meta$`Sequencing approach`[which(phage.meta$`Sequence ID a` == x)]))
sapply(ARGlist,class)

ARGlist$Family[ARGlist$Family == "n.a."] <- "Unidentified"


# total number of ARG in each bar  -----
totalARGnum <- ARGlist %>%
  dplyr::group_by(Family, `Sequencing approach`) %>%
  summarise(numTotalARG = n()) %>% 
  as.data.frame()



# fraction of ARG types in each genome -----
tmp <- ARGlist %>%
  group_by(`Sequence ID`, `Sequencing approach`, Family, `Antibiotic type`) %>%
  summarise(numORF = n()) %>%
  as.data.frame() %>%
  group_by(Family, `Sequencing approach`, `Sequence ID`) %>%
  mutate(perc = numORF/sum(numORF)) %>% # the percentage of each ARG type in each genome
  as.data.frame()
tmp %>% group_by( `Sequence ID`) %>% summarise(s=sum(perc)) # check results 



# complete the 0 values
tmp2 <- merge(tmp %>% 
                reshape2::dcast(`Sequence ID` ~ `Antibiotic type`, value.var = "perc") %>%
                reshape2::melt(variable.name="Antibiotic type"),
              tmp %>% select(`Sequence ID`, `Sequencing approach`, Family) %>% unique,
              by=c("Sequence ID"), all = T) 
tmp2$value[is.na(tmp2$value)] <- 0 # value是percentage


argPerc.perGenome_phage <- merge(tmp2, totalARGnum, by=c("Sequencing approach","Family") ) 



#  fraction of ARG types each bar
plotDat.phage <- merge(tmp2 %>%
                         group_by(Family, `Sequencing approach`, `Antibiotic type`) %>%
                         summarise(avg.perc= mean(value)) %>%
                         as.data.frame(),
                       totalARGnum, 
                       by=c("Family","Sequencing approach")) %>%
  mutate(bar = paste(`Sequencing approach`, Family, sep = "|")) 

plotDat.phage <- plotDat.phage %>% filter(bar %in% kept.bars)
plotDat.phage %>% group_by(Family, `Sequencing approach`) %>% summarise(s=sum(avg.perc)) # check results 



# add overall to the plotDat 
tmp.overall<-  
  merge(
    tmp2 %>%
      group_by( `Sequencing approach`, `Antibiotic type`) %>%
      summarise(avg.perc= mean(value)) %>%
      as.data.frame() %>%
      mutate(Family="Overall"),
    totalARGnum %>%
      group_by( `Sequencing approach`) %>%
      summarise(numTotalARG= sum(numTotalARG)) %>%
      as.data.frame() %>%
      mutate(Family="Overall"),
    by=c("Family","Sequencing approach")
  )

plotDat.phage <- bind_rows(plotDat.phage, tmp.overall)


# colors ----
plotDat.phage$`Antibiotic type` <- as.character(plotDat.phage$`Antibiotic type`)
ARGtype.lvl <- unique(plotDat.phage$`Antibiotic type`)[order(unique(plotDat.phage$`Antibiotic type`))]
colors <- sapply(ARGtype.lvl,
                 function(x){
                   colors.ARGtype$color[colors.ARGtype$ARGtype == x] }) 





# 按IG和MAG分别画图 ----

for(tp in c("Isolate","Metagenomic")){
  # tp="Metagenomic"  
  
  plotDat.phage_sub <- plotDat.phage %>% filter(`Sequencing approach` == tp)
  fam.lvls <- eval(parse(text = paste0("fam.lvl_", tp)))
  
  plotDat.phage_sub$Family <- factor(plotDat.phage_sub$Family,
                                     levels = fam.lvls)
  
  # 不要太多的colors:
  
  colors.sub <- colors[unique((plotDat.phage_sub %>%
                                 filter(avg.perc > 0))$`Antibiotic type`)]
  
  P <- ggplot(plotDat.phage_sub) +
    #geom_alluvium(aes(fill = `Antibiotic type`),alpha = 0.6, width = 0.6) +
    #geom_stratum(aes(fill = `Antibiotic type`),width = 0.6, color="white", size=0.15) +
    geom_col(aes(x=Family, y=avg.perc, fill=`Antibiotic type`), width = 0.6, color="white", size=0.15 ) + 
    theme_bw() + 
    theme( # axis.text.x = element_text(angle = 90),
      panel.grid = element_blank()) +
    scale_fill_manual(values = colors.sub)   +  
    # theme(legend.position = "none") +
    scale_y_continuous(labels = scales::percent) 
  P
  assign(paste0("P.phage_",tp), P)
  
}



ggarrange(P.phage_Isolate + theme(legend.position = "none", axis.text.x = element_text(angle = 90)),
          P.phage_Metagenomic + theme(legend.position = "none", axis.text.x = element_text(angle = 90)), 
          widths = c(0.59, 0.41))

pdf("Fig2D.ARGtypeComposition_phage.pdf",  width = 3.3, height = 1.9)
ggarrange(P.phage_Isolate + theme(legend.position = "none"),
          P.phage_Metagenomic + theme(legend.position = "none"), widths = c(0.59, 0.41))
dev.off()
