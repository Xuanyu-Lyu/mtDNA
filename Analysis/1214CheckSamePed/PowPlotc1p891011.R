library(ggplot2)

powerCal <- function(lamda,df){
      1- pchisq(qchisq(1-.05, df), df, lamda)
      
}
my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")
n = 2000
#relativePath = "C:/Users/lxy75/OneDrive/Documents/R-Project/mtDNA/Data/1130Prelim"
#relativePath_save = "C:/Users/lxy75/OneDrive/Documents/R-Project/mtDNA/Analysis/1130Prelim"

relativePath = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/1214CheckSamePed"
relativePath_save = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Analysis/1214CheckSamePed"

### create a new enviroment for one condition
env_c1p8 <- new.env()
load(paste0(relativePath,"/","c1p8/modelSmr.Rdata"), envir = env_c1p8)


meanDiffLL_mtam_c1p8 <- env_c1p8$smr2$Minus2LogLikelihood - env_c1p8$smr1$Minus2LogLikelihood
PN <- ncol(env_c1p8[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c1p8[["smr1"]][["dataSummary"]])
lamdaUnit_c1p8 <- meanDiffLL_mtam_c1p8/PN
SSize_c1p8 <- 1: n
LamdaVec_c1p8 <- lamdaUnit_c1p8*SSize_c1p8
powVec_c1p8 <- as.numeric(lapply(LamdaVec_c1p8,powerCal, df = 1))

df_c1p8 <- data.frame(Nped = SSize_c1p8, 
                      power = powVec_c1p8, 
                      Combination = rep("power1", n))


### create a new enviroment for one condition
env_c1p9 <- new.env()
load(paste0(relativePath,"/","c1p9/modelSmr.Rdata"), envir = env_c1p9)


meanDiffLL_mtam_c1p9 <- env_c1p9$smr2$Minus2LogLikelihood - env_c1p9$smr1$Minus2LogLikelihood
PN <- ncol(env_c1p9[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c1p9[["smr1"]][["dataSummary"]])
lamdaUnit_c1p9 <- meanDiffLL_mtam_c1p9/PN
SSize_c1p9 <- 1: n
LamdaVec_c1p9 <- lamdaUnit_c1p9*SSize_c1p9
powVec_c1p9 <- as.numeric(lapply(LamdaVec_c1p9,powerCal, df = 1))

df_c1p9 <- data.frame(Nped = SSize_c1p9, 
                      power = powVec_c1p9, 
                      Combination = rep("power2", n))


### create a new enviroment for one condition
env_c1p10 <- new.env()
load(paste0(relativePath,"/","c1p10/modelSmr.Rdata"), envir = env_c1p10)


meanDiffLL_mtam_c1p10 <- env_c1p10$smr2$Minus2LogLikelihood - env_c1p10$smr1$Minus2LogLikelihood
PN <- ncol(env_c1p10[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c1p10[["smr1"]][["dataSummary"]])
lamdaUnit_c1p10 <- meanDiffLL_mtam_c1p10/PN
SSize_c1p10 <- 1: n
LamdaVec_c1p10 <- lamdaUnit_c1p10*SSize_c1p10
powVec_c1p10 <- as.numeric(lapply(LamdaVec_c1p10,powerCal, df = 1))

df_c1p10 <- data.frame(Nped = SSize_c1p10, 
                      power = powVec_c1p10, 
                      Combination = rep("power3", n))

### create a new enviroment for one condition
env_c1p11 <- new.env()
load(paste0(relativePath,"/","c1p11/modelSmr.Rdata"), envir = env_c1p11)


meanDiffLL_mtam_c1p11 <- env_c1p11$smr2$Minus2LogLikelihood - env_c1p11$smr1$Minus2LogLikelihood
PN <- ncol(env_c1p11[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c1p11[["smr1"]][["dataSummary"]])
lamdaUnit_c1p11 <- meanDiffLL_mtam_c1p11/PN
SSize_c1p11 <- 1: n
LamdaVec_c1p11 <- lamdaUnit_c1p11*SSize_c1p11
powVec_c1p11 <- as.numeric(lapply(LamdaVec_c1p11,powerCal, df = 1))

df_c1p11 <- data.frame(Nped = SSize_c1p11, 
                      power = powVec_c1p11, 
                      Combination = rep("power4", n))

# ### create a new enviroment for one condition
# env_c1p5 <- new.env()
# load(paste0(relativePath,"/","c1p5/modelSmr.Rdata"), envir = env_c1p5)


# meanDiffLL_mtam_c1p5 <- env_c1p5$smr2$Minus2LogLikelihood - env_c1p5$smr1$Minus2LogLikelihood
# PN <- ncol(env_c1p5[["smr1"]][["dataSummary"]][["fam1"]])*length(env_c1p5[["smr1"]][["dataSummary"]])
# lamdaUnit_c1p5 <- meanDiffLL_mtam_c1p5/PN
# SSize_c1p5 <- 1: n
# LamdaVec_c1p5 <- lamdaUnit_c1p5*SSize_c1p5
# powVec_c1p5 <- as.numeric(lapply(LamdaVec_c1p5,powerCal, df = 1))

# df_c1p5 <- data.frame(Nped = SSize_c1p5, 
#                       power = powVec_c1p5, 
#                       Combination = rep("power5", n))

### create a data frame for graphs
df_c1 <- rbind(df_c1p8, df_c1p9, df_c1p10,df_c1p11)
df_c1$Combination <- as.factor(df_c1$Combination)


g1 <-ggplot(data = df_c1)+ geom_line(mapping = aes(x = Nped, y = power, color= Combination), size = 1.5) +
      scale_color_manual(values=my_palette[1:4],
                         name="Pedigree Structures",
                         breaks=c("power1", "power2", "power3", "power4"),
                         labels=c("k = 4, G = 5, m = NA",
                                  "k = 4, G = 5, m = NA", 
                                  "k = 4, G = 5, m = NA",
                                  "k = 4, G = 5, m = NA")
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

ggsave( paste0(relativePath_save,"/graph_c1p891011.png"), g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)

#rm(list = ls())

