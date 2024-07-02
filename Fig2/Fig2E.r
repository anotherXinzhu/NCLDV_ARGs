
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(scales)
library(RColorBrewer )


# colors =====================================
library(paletteer) 

colors.resMech <- 
  cbind.data.frame(
    resMech = c( "Antibiotic inactivation","Antibiotic efflux","Antibiotic target alteration",
                 #"Antibiotic target protection","Antibiotic target replacement",
                 "Antibiotic efflux;Antibiotic target alteration","unknown"),
    color = c(paletteer_d("ggthemes::Red_Blue_Brown")[1:3],
              "#9D7660FF","gray"),
    stringsAsFactors=F
  )


# the family levels -----
load("../Fig1/Fig1AD.plotSourceDat.RData") # for the meta data
load("Fig2A.plots_n.sourceDat_ncldv.RData") # this is for the levels of kept bars 
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

i.correct <- which(ARGlist$`Resistance mechanism` %in% c("Antibiotic target alteration", "Antibiotic target protection","Antibiotic target replacement"))
ARGlist$`Resistance mechanism`[i.correct] <- "Antibiotic target alteration"


# total number of ARG in each bar  -----
totalARGnum <- ARGlist %>%
  dplyr::group_by(Family, `Sequencing approach`) %>%
  summarise(numTotalARG = n()) %>% 
  as.data.frame()


# fraction of ARG types in each genome -----
tmp <- ARGlist %>%
  group_by(`Sequence ID`, `Sequencing approach`, Family, `Resistance mechanism`) %>%
  summarise(numORF = n()) %>%
  as.data.frame() %>%
  group_by(Family, `Sequencing approach`, `Sequence ID`) %>%
  mutate(perc = numORF/sum(numORF)) %>% # the percentage of each ARG type in each genome
  as.data.frame()
tmp %>% group_by(`Sequence ID`) %>% summarise(s=sum(perc)) # check results 

# complete the 0 values
tmp2 <- merge(tmp %>% 
                reshape2::dcast(`Sequence ID` ~ `Resistance mechanism`, value.var = "perc") %>%
                reshape2::melt(variable.name="Resistance mechanism"),
              tmp %>% select(`Sequence ID`, `Sequencing approach`, Family) %>% unique,
              by=c("Sequence ID"), all = T) 
tmp2$value[is.na(tmp2$value)] <- 0


argPerc.perGenome_ncldv <- merge(tmp2, totalARGnum, by=c("Sequencing approach","Family") ) # to calculate statistical significance



plotDat.ncldv <- merge(tmp2 %>%
                         group_by(Family, `Sequencing approach`, `Resistance mechanism`) %>%
                         summarise(avg.perc= mean(value)) %>%
                         as.data.frame(),
                       totalARGnum, 
                       by=c("Family","Sequencing approach")) %>%
  mutate(bar = paste(`Sequencing approach`, Family, sep = "|")) 

plotDat.ncldv <- plotDat.ncldv %>% filter(bar %in% kept.bars)
plotDat.ncldv %>% group_by(Family, `Sequencing approach`) %>% summarise(s=sum(avg.perc)) # check results 


# add overall to the plotDat 
tmp.overall <-  
  merge(
    tmp2 %>%
      group_by( `Sequencing approach`, `Resistance mechanism`) %>%
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
plotDat.ncldv$`Resistance mechanism` <- as.character(plotDat.ncldv$`Resistance mechanism`)
ARGtype.lvl <- unique(plotDat.ncldv$`Resistance mechanism`)[order(unique(plotDat.ncldv$`Resistance mechanism`))]
colors <- sapply(ARGtype.lvl,
                 function(x){
                   colors.resMech$color[colors.resMech$resMech == x] }) 


# plotting ------

for(tp in c("Isolate","Metagenomic")){
  # tp = "Isolate"  
  
  plotDat.ncldv_sub <- plotDat.ncldv %>% filter(`Sequencing approach` == tp)
  fam.lvls <- eval(parse(text = paste0("fam.lvl_", tp)))
  
  plotDat.ncldv_sub$Family <- factor(plotDat.ncldv_sub$Family,
                                     levels = fam.lvls)
  
  colors.sub <- colors[unique((plotDat.ncldv_sub %>%
                                 filter(avg.perc > 0))$`Resistance mechanism`)]
  
  P <- ggplot(plotDat.ncldv_sub) +
    geom_col(aes(x=Family, y=avg.perc, fill=`Resistance mechanism`), width = 0.6, color="black", size=0.15 ) + 
    theme_bw() + 
    theme( # axis.text.x = element_text(angle = 90),
      panel.grid = element_blank()) +
    scale_fill_manual(values = colors.sub)   +  
    scale_y_continuous(labels = scales::percent) 
  P
  assign(paste0("P.ncldv_",tp), P)
  
}



pdf("Fig2E.resMechanisms_ncldvs.pdf",   width = 5, height = 1.9)
ggarrange(P.ncldv_Isolate + theme(legend.position = "none"),
          P.ncldv_Metagenomic + theme(legend.position = "none"), widths = c(0.56, 0.44))
dev.off()

