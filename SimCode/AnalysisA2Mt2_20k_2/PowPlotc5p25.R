library(ggplot2)

powerCal <- function(lamda,df){
    1- pchisq(qchisq(1-.05, df), df, lamda)
    
}
my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")


### create a new enviroment for one condition
env_c5p2 <- new.env()
load("~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/c5p2/modelSmr.Rdata", envir = env_c5p2)


meanDiffLL_mtam_c5p2 <- env_c5p2$smr2$Minus2LogLikelihood - env_c5p2$smr1$Minus2LogLikelihood
PN <- ncol(env_c5p2[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c5p2[["smr1"]][["dataSummary"]])
lamdaUnit_c5p2 <- meanDiffLL_mtam_c5p2/PN
SSize_c5p2 <- 1: 2500
LamdaVec_c5p2 <- lamdaUnit_c5p2*SSize_c5p2
powVec_c5p2 <- as.numeric(lapply(LamdaVec_c5p2,powerCal, df = 1))

df_c5p2 <- data.frame(Nped = SSize_c5p2, 
                      power = powVec_c5p2, 
                      Combination = rep("power1", 2500))


### create a new enviroment for one condition
env_c5p5 <- new.env()
load("~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/c5p5/modelSmr.Rdata", envir = env_c5p5)


meanDiffLL_mtam_c5p5 <- env_c5p5$smr2$Minus2LogLikelihood - env_c5p5$smr1$Minus2LogLikelihood
PN <- ncol(env_c5p5[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c5p5[["smr1"]][["dataSummary"]])
lamdaUnit_c5p5 <- meanDiffLL_mtam_c5p5/PN
SSize_c5p5 <- 1: 2500
LamdaVec_c5p5 <- lamdaUnit_c5p5*SSize_c5p5
powVec_c5p5 <- as.numeric(lapply(LamdaVec_c5p5,powerCal, df = 1))

df_c5p5 <- data.frame(Nped = SSize_c5p5, 
                      power = powVec_c5p5, 
                      Combination = rep("power2", 2500))




### create a data frame for graphs
df_c5 <- rbind(df_c5p2, df_c5p5)
df_c5$Combination <- as.factor(df_c5$Combination)


g1 <-ggplot(data = df_c5)+ geom_line(mapping = aes(x = Nped, y = power, color= Combination), size = 1.5) +
    scale_color_manual(values=my_palette[1:2],
                       name="Pedigree Structures",
                       breaks=c("power1", "power2" ),
                       labels=c("k = 2, G = 4, m = 17",
                                "k = 2, G = 3, m = 9")
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

ggsave( "~/R-Project/mtDNA_mt2/ResultA2Mt2_20k/AnalysisA2Mt2_20k/graph_c5p25.png",g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)

#rm(list = ls())
