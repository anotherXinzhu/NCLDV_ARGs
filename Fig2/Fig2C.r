
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
    ARGtype = c( "aminocoumarin","phenicol","multidrug","polymyxin","tetracycline","trimethoprim","mupirocin", 
                 "sulfonamide","peptide","MLS","rifamycin","disinfecting agents and antiseptics","quinolone","nitroimidazole","pleuromutilin","Others"),
    color = c(getPalette(colourCount), "gray"),
    stringsAsFactors=F
  )

# the family levels -----
load("../Fig1/Fig1AD.plotSourceDat.RData") # for the meta data
load("Fig2A.plots_n.sourceDat_ncldv.RData") # this is for the levels of kept bars and family levels
fam.lvl_Isolate <- (avg.plotdat.ncldv_Isolate %>% arrange(desc(avg.type)) %>% arrange(group.type))$Family
fam.lvl_Metagenomic <- (avg.plotdat.ncldv_Metagenomic %>% arrange(desc(avg.type)) %>% arrange(group.type))$Family
remove(ARGdiver_df)
remove(list = ls()[grepl("^P", ls())])
remove(list = ls()[grepl("^testDat", ls())])
remove(list = ls()[grepl("^avg.plotdat.ncldv", ls())])


# the ARG list ----
library(readxl)
ARGlist <- read_excel("../_data/TableS.xlsx", sheet = 4, skip= 2)
as.character(ARGlist[1,4:6])
colnames(ARGlist)[4:6] <- as.character(ARGlist[1,4:6])
ARGlist <- ARGlist[2:750,]

ARGlist <- ARGlist %>%
  mutate(Family = sapply(`Sequence ID`, function(x) ncldv.meta$Family[which(ncldv.meta$`Sequence ID a`== x)])) %>%
  mutate(`Sequencing approach` = sapply( `Sequence ID`, function(x)ncldv.meta$`Sequencing approach`[which(ncldv.meta$`Sequence ID a` == x)]))
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
tmp %>% group_by(`Sequence ID`) %>% summarise(s=sum(perc)) # check results 



# complete the 0 values
tmp2 <- merge(tmp %>% 
                reshape2::dcast(`Sequence ID` ~ `Antibiotic type`, value.var = "perc") %>%
                reshape2::melt(variable.name="Antibiotic type"),
              tmp %>% select(`Sequence ID`, `Sequencing approach`, Family) %>% unique,
              by=c("Sequence ID"), all = T) 
tmp2$value[is.na(tmp2$value)] <- 0 # value is percentage


argPerc.perGenome_ncldv <- merge(tmp2, totalARGnum, by=c("Sequencing approach","Family") ) # to calculate statistical significance


#  fraction of ARG types each bar
plotDat.ncldv <- merge(tmp2 %>%
                         group_by(Family, `Sequencing approach`, `Antibiotic type`) %>%
                         summarise(avg.perc= mean(value)) %>%
                         as.data.frame(),
                       totalARGnum, 
                       by=c("Family","Sequencing approach")) %>%
  mutate(bar = paste(`Sequencing approach`, Family, sep = "|")) 

plotDat.ncldv <- plotDat.ncldv %>% filter(bar %in% kept.bars)
plotDat.ncldv %>% group_by(Family, `Sequencing approach`) %>% summarise(s=sum(avg.perc)) # check results 



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

plotDat.ncldv <- bind_rows(plotDat.ncldv, tmp.overall)


# colors ----
plotDat.ncldv$`Antibiotic type` <- as.character(plotDat.ncldv$`Antibiotic type`)
ARGtype.lvl <- unique(plotDat.ncldv$`Antibiotic type`)[order(unique(plotDat.ncldv$`Antibiotic type`))]
colors <- sapply(ARGtype.lvl,
                 function(x){
                   colors.ARGtype$color[colors.ARGtype$ARGtype == x] }) 


# plotting ----

for(tp in c("Isolate","Metagenomic")){
  # tp="Isolate"  
  
  plotDat.ncldv_sub <- plotDat.ncldv %>% filter(`Sequencing approach` == tp)
  fam.lvls <- eval(parse(text = paste0("fam.lvl_", tp)))
  
  plotDat.ncldv_sub$Family <- factor(plotDat.ncldv_sub$Family,
                                     levels = fam.lvls)
  
  colors.sub <- colors[unique((plotDat.ncldv_sub %>%
                                 filter(avg.perc > 0))$`Antibiotic type`)]
  
  P <- ggplot(plotDat.ncldv_sub) +
    geom_col(aes(x=Family, y=avg.perc, fill=`Antibiotic type`), width = 0.6, color="white", size=0.15 ) + 
    theme_bw() + 
    theme( # axis.text.x = element_text(angle = 90),
      panel.grid = element_blank()) +
    scale_fill_manual(values = colors.sub)   +  
    scale_y_continuous(labels = scales::percent) 
  P
  assign(paste0("P.ncldv_",tp), P)
  
}

ggarrange(P.ncldv_Isolate + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) ,
          P.ncldv_Metagenomic + theme(legend.position = "none", axis.text.x = element_text(angle = 90)), 
          widths = c(0.56, 0.44))

pdf("Fig2C.ARGtypeComposition_ncldv.pdf",  width = 5, height = 1.9)
ggarrange(P.ncldv_Isolate + theme(legend.position = "none") ,
          P.ncldv_Metagenomic + theme(legend.position = "none"), 
          widths = c(0.56, 0.44))
dev.off()
