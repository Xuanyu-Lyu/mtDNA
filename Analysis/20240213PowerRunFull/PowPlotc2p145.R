
library(ggplot2)

powerCal <- function(lamda,df){
      1- pchisq(qchisq(1-.05, df), df, lamda)
      
}
my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")
n = 10000
#relativePath = "C:/Users/lxy75/OneDrive/Documents/R-Project/mtDNA/Data/1130Prelim"
#relativePath_save = "C:/Users/lxy75/OneDrive/Documents/R-Project/mtDNA/Analysis/1130Prelim"

relativePath = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/20240213PowerRunFull"
relativePath_save = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Analysis/20240213PowerRunFull"

### create a new enviroment for one condition
env_c2p1 <- new.env()
load(paste0(relativePath,"/","c2p1/modelSmr.Rdata"), envir = env_c2p1)


meanDiffLL_mtam_c2p1 <- env_c2p1$smr2$Minus2LogLikelihood - env_c2p1$smr1$Minus2LogLikelihood
PN <- ncol(env_c2p1[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p1[["smr1"]][["dataSummary"]])
lamdaUnit_c2p1 <- meanDiffLL_mtam_c2p1/PN
SSize_c2p1 <- 1: n
LamdaVec_c2p1 <- lamdaUnit_c2p1*SSize_c2p1
powVec_c2p1 <- as.numeric(lapply(LamdaVec_c2p1,powerCal, df = 2))

df_c2p1 <- data.frame(Nped = SSize_c2p1, 
                      power = powVec_c2p1, 
                      Combination = rep("power1", n))


### create a new enviroment for one condition
env_c2p4 <- new.env()
load(paste0(relativePath,"/","c2p4/modelSmr.Rdata"), envir = env_c2p4)


meanDiffLL_mtam_c2p4 <- env_c2p4$smr2$Minus2LogLikelihood - env_c2p4$smr1$Minus2LogLikelihood
PN <- ncol(env_c2p4[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p4[["smr1"]][["dataSummary"]])
lamdaUnit_c2p4 <- meanDiffLL_mtam_c2p4/PN
SSize_c2p4 <- 1: n
LamdaVec_c2p4 <- lamdaUnit_c2p4*SSize_c2p4
powVec_c2p4 <- as.numeric(lapply(LamdaVec_c2p4,powerCal, df = 2))

df_c2p4 <- data.frame(Nped = SSize_c2p4, 
                      power = powVec_c2p4, 
                      Combination = rep("power2", n))


### create a new enviroment for one condition
env_c2p5 <- new.env()
load(paste0(relativePath,"/","c2p5/modelSmr.Rdata"), envir = env_c2p5)


meanDiffLL_mtam_c2p5 <- env_c2p5$smr2$Minus2LogLikelihood - env_c2p5$smr1$Minus2LogLikelihood
PN <- ncol(env_c2p5[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p5[["smr1"]][["dataSummary"]])
lamdaUnit_c2p5 <- meanDiffLL_mtam_c2p5/PN
SSize_c2p5 <- 1: n
LamdaVec_c2p5 <- lamdaUnit_c2p5*SSize_c2p5
powVec_c2p5 <- as.numeric(lapply(LamdaVec_c2p5,powerCal, df = 2))

df_c2p5 <- data.frame(Nped = SSize_c2p5, 
                      power = powVec_c2p5, 
                      Combination = rep("power3", n))

# ### create a new enviroment for one condition
# env_c2p4 <- new.env()
# load(paste0(relativePath,"/","c2p4/modelSmr.Rdata"), envir = env_c2p4)


# meanDiffLL_mtam_c2p4 <- env_c2p4$smr2$Minus2LogLikelihood - env_c2p4$smr1$Minus2LogLikelihood
# PN <- ncol(env_c2p4[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p4[["smr1"]][["dataSummary"]])
# lamdaUnit_c2p4 <- meanDiffLL_mtam_c2p4/PN
# SSize_c2p4 <- 1: n
# LamdaVec_c2p4 <- lamdaUnit_c2p4*SSize_c2p4
# powVec_c2p4 <- as.numeric(lapply(LamdaVec_c2p4,powerCal, df = 1))

# df_c2p4 <- data.frame(Nped = SSize_c2p4, 
#                       power = powVec_c2p4, 
#                       Combination = rep("power4", n))

# ### create a new enviroment for one condition
# env_c2p5 <- new.env()
# load(paste0(relativePath,"/","c2p5/modelSmr.Rdata"), envir = env_c2p5)


# meanDiffLL_mtam_c2p5 <- env_c2p5$smr2$Minus2LogLikelihood - env_c2p5$smr1$Minus2LogLikelihood
# PN <- ncol(env_c2p5[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c2p5[["smr1"]][["dataSummary"]])
# lamdaUnit_c2p5 <- meanDiffLL_mtam_c2p5/PN
# SSize_c2p5 <- 1: n
# LamdaVec_c2p5 <- lamdaUnit_c2p5*SSize_c2p5
# powVec_c2p5 <- as.numeric(lapply(LamdaVec_c2p5,powerCal, df = 1))

# df_c2p5 <- data.frame(Nped = SSize_c2p5, 
#                       power = powVec_c2p5, 
#                       Combination = rep("power5", n))

### create a data frame for graphs
df_c2 <- rbind(df_c2p1, df_c2p4, df_c2p5
               #,df_c2p4,df_c2p5
               )
df_c2$Combination <- as.factor(df_c2$Combination)


g1 <-ggplot(data = df_c2)+ geom_line(mapping = aes(x = Nped, y = power, color= Combination), size = 1.5) +
      scale_color_manual(values=my_palette[1:3],
                         name="Pedigree Structures",
                         breaks=c("power1", "power2", "power3" 
                                  #,"power4","power5"
                                  ),
                         labels=c("k = 3, G = 4, m = 33",
                                  "k = 3, G = 6, m = 388",
                                  "k = 3, G = 8, m = 657")
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
      ylab("Power (mt\u00B2)")+
      geom_hline(yintercept = .8, linetype = 5, size = .8, color = "grey")
# +
# annotate(geom = "text",x = 0.62, y =.92, label = "a\u00B2 = .6", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.73, y =.6, label = "a\u00B2 = .4", family="Calibri", color = "gray40",size = 3)+
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .2", family="Calibri", color = "gray40",size = 3)+ 
# annotate(geom = "text",x = 0.9, y =.6, label = "a\u00B2 = .4, d\u00B2 = .1", family="Calibri", color = "gray40",size = 3)
g1

ggsave( paste0(relativePath_save,"/graph_c2p145.png"), g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)

#rm(list = ls())
