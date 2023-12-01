library(ggplot2)

powerCal <- function(lamda,df){
    1- pchisq(qchisq(1-.05, df), df, lamda)
    
}
my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")


### create a new enviroment for one condition
env_c1p1 <- new.env()
load("~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/c1p1/modelSmr.Rdata", envir = env_c1p1)


meanDiffLL_mtam_c1p1 <- env_c1p1$smr2$Minus2LogLikelihood - env_c1p1$smr1$Minus2LogLikelihood
PN <- ncol(env_c1p1[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c1p1[["smr1"]][["dataSummary"]])
lamdaUnit_c1p1 <- meanDiffLL_mtam_c1p1/PN
SSize_c1p1 <- 1: 10000
LamdaVec_c1p1 <- lamdaUnit_c1p1*SSize_c1p1
powVec_c1p1 <- as.numeric(lapply(LamdaVec_c1p1,powerCal, df = 1))

df_c1p1 <- data.frame(Nped = SSize_c1p1, 
                      power = powVec_c1p1, 
                      Combination = rep("power1", 10000))


### create a new enviroment for one condition
env_c1p2 <- new.env()
load("~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/c1p2/modelSmr.Rdata", envir = env_c1p2)


meanDiffLL_mtam_c1p2 <- env_c1p2$smr2$Minus2LogLikelihood - env_c1p2$smr1$Minus2LogLikelihood
PN <- ncol(env_c1p2[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c1p2[["smr1"]][["dataSummary"]])
lamdaUnit_c1p2 <- meanDiffLL_mtam_c1p2/PN
SSize_c1p2 <- 1: 10000
LamdaVec_c1p2 <- lamdaUnit_c1p2*SSize_c1p2
powVec_c1p2 <- as.numeric(lapply(LamdaVec_c1p2,powerCal, df = 1))

df_c1p2 <- data.frame(Nped = SSize_c1p2, 
                      power = powVec_c1p2, 
                      Combination = rep("power2", 10000))


### create a new enviroment for one condition
env_c1p3 <- new.env()
load("~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/c1p3/modelSmr.Rdata", envir = env_c1p3)


meanDiffLL_mtam_c1p3 <- env_c1p3$smr2$Minus2LogLikelihood - env_c1p3$smr1$Minus2LogLikelihood
PN <- ncol(env_c1p3[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c1p3[["smr1"]][["dataSummary"]])
lamdaUnit_c1p3 <- meanDiffLL_mtam_c1p3/PN
SSize_c1p3 <- 1: 10000
LamdaVec_c1p3 <- lamdaUnit_c1p3*SSize_c1p3
powVec_c1p3 <- as.numeric(lapply(LamdaVec_c1p3,powerCal, df = 1))

df_c1p3 <- data.frame(Nped = SSize_c1p3, 
                      power = powVec_c1p3, 
                      Combination = rep("power3", 10000))

### create a new enviroment for one condition
env_c1p4 <- new.env()
load("~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/c1p4/modelSmr.Rdata", envir = env_c1p4)


meanDiffLL_mtam_c1p4 <- env_c1p4$smr2$Minus2LogLikelihood - env_c1p4$smr1$Minus2LogLikelihood
PN <- ncol(env_c1p4[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c1p4[["smr1"]][["dataSummary"]])
lamdaUnit_c1p4 <- meanDiffLL_mtam_c1p4/PN
SSize_c1p4 <- 1: 10000
LamdaVec_c1p4 <- lamdaUnit_c1p4*SSize_c1p4
powVec_c1p4 <- as.numeric(lapply(LamdaVec_c1p4,powerCal, df = 1))

df_c1p4 <- data.frame(Nped = SSize_c1p4, 
                      power = powVec_c1p4, 
                      Combination = rep("power4", 10000))

### create a data frame for graphs
df_c1 <- rbind(df_c1p1, df_c1p2, df_c1p3,df_c1p4)
df_c1$Combination <- as.factor(df_c1$Combination)


g1 <-ggplot(data = df_c1)+ geom_line(mapping = aes(x = Nped, y = power, color= Combination), size = 1.5) +
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

ggsave( "~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/AnalysisA2Mt2_20k/graph_c1p1234.png",g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)

#rm(list = ls())
