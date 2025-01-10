library(ggplot2)

powerCal <- function(lamda,df){
    1- pchisq(qchisq(1-.05, df), df, lamda)
    
}
my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")


### create a new enviroment for one condition
env_c4p1 <- new.env()
load("~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/c4p1/modelSmr.Rdata", envir = env_c4p1)


meanDiffLL_mtam_c4p1 <- env_c4p1$smr2$Minus2LogLikelihood - env_c4p1$smr1$Minus2LogLikelihood
PN <- ncol(env_c4p1[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c4p1[["smr1"]][["dataSummary"]])
lamdaUnit_c4p1 <- meanDiffLL_mtam_c4p1/PN
SSize_c4p1 <- 1: 10000
LamdaVec_c4p1 <- lamdaUnit_c4p1*SSize_c4p1
powVec_c4p1 <- as.numeric(lapply(LamdaVec_c4p1,powerCal, df = 1))

df_c4p1 <- data.frame(Nped = SSize_c4p1, 
                      power = powVec_c4p1, 
                      Combination = rep("power1", 10000))


### create a new enviroment for one condition
env_c4p2 <- new.env()
load("~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/c4p2/modelSmr.Rdata", envir = env_c4p2)


meanDiffLL_mtam_c4p2 <- env_c4p2$smr2$Minus2LogLikelihood - env_c4p2$smr1$Minus2LogLikelihood
PN <- ncol(env_c4p2[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c4p2[["smr1"]][["dataSummary"]])
lamdaUnit_c4p2 <- meanDiffLL_mtam_c4p2/PN
SSize_c4p2 <- 1: 10000
LamdaVec_c4p2 <- lamdaUnit_c4p2*SSize_c4p2
powVec_c4p2 <- as.numeric(lapply(LamdaVec_c4p2,powerCal, df = 1))

df_c4p2 <- data.frame(Nped = SSize_c4p2, 
                      power = powVec_c4p2, 
                      Combination = rep("power2", 10000))


### create a new enviroment for one condition
env_c4p3 <- new.env()
load("~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/c4p3/modelSmr.Rdata", envir = env_c4p3)


meanDiffLL_mtam_c4p3 <- env_c4p3$smr2$Minus2LogLikelihood - env_c4p3$smr1$Minus2LogLikelihood
PN <- ncol(env_c4p3[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c4p3[["smr1"]][["dataSummary"]])
lamdaUnit_c4p3 <- meanDiffLL_mtam_c4p3/PN
SSize_c4p3 <- 1: 10000
LamdaVec_c4p3 <- lamdaUnit_c4p3*SSize_c4p3
powVec_c4p3 <- as.numeric(lapply(LamdaVec_c4p3,powerCal, df = 1))

df_c4p3 <- data.frame(Nped = SSize_c4p3, 
                      power = powVec_c4p3, 
                      Combination = rep("power3", 10000))

### create a new enviroment for one condition
env_c4p4 <- new.env()
load("~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/c4p4/modelSmr.Rdata", envir = env_c4p4)


meanDiffLL_mtam_c4p4 <- env_c4p4$smr2$Minus2LogLikelihood - env_c4p4$smr1$Minus2LogLikelihood
PN <- ncol(env_c4p4[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c4p4[["smr1"]][["dataSummary"]])
lamdaUnit_c4p4 <- meanDiffLL_mtam_c4p4/PN
SSize_c4p4 <- 1: 10000
LamdaVec_c4p4 <- lamdaUnit_c4p4*SSize_c4p4
powVec_c4p4 <- as.numeric(lapply(LamdaVec_c4p4,powerCal, df = 1))

df_c4p4 <- data.frame(Nped = SSize_c4p4, 
                      power = powVec_c4p4, 
                      Combination = rep("power4", 10000))

### create a data frame for graphs
df_c4 <- rbind(df_c4p1, df_c4p2, df_c4p3,df_c4p4)
df_c4$Combination <- as.factor(df_c4$Combination)


g1 <-ggplot(data = df_c4)+ geom_line(mapping = aes(x = Nped, y = power, color= Combination), size = 1.5) +
    scale_color_manual(values=my_palette[1:4],
                       name="Pedigree Structures",
                       breaks=c("power1", "power2", "power3", "power4"),
                       labels=c("k = 2, G = 4, m = 17",
                                "k = 3, G = 4, m = 36", 
                                "k = 4, G = 4, m = 66",
                                "k = 8, G = 4, m = 388")
    )+
    theme(panel.background = element_rect(fill = "transparent"),
          panel.grid = element_line(color = "transparent"),
          axis.line = element_line(size = 1, colour = "black"),
          #axis.line.y = element_blank(),
          axis.text = element_text( color = "black"),
          #axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          text=element_text( family="Calibri",  size = 12),
          legend.spacing = unit(-17,'pt'),
          legend.margin = margin(t=0,b=0,unit='pt'),
          legend.background = element_rect(),
          legend.position=c(.8,.2))+
    xlab("N of Individuals")+
    scale_y_continuous(n.breaks = 6)+
    ylab("Power (mt\u00B2+j\u00B2)")+
    geom_hline(yintercept = .8, linetype = 5, size = .8, color = "grey")
# +
# annotate(geom = "text",x = 0.62, y =.92, label = "a\u00B2 = .6", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.73, y =.6, label = "a\u00B2 = .4", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .2", family="Calibri", color = "gray40",size = 3)+ 
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .4, d\u00B2 = .1", family="Calibri", color = "gray40",size = 3)
g1

ggsave( "~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/AnalysisA2Mt2_20k/graph_c4p1234.png",g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)

#rm(list = ls())
