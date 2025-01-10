
library(ggplot2)

powerCal <- function(lamda,df){
      1- pchisq(qchisq(1-.05, df), df, lamda)
      
}
my_palette <- c( "#E41A1C", "#332288", "#E69F00", "#DDCC77", "#377EB8",  "#4DAF4A", "#117A65", "#56B4E9", "#A6CE39", "#A9A9A9","#88CCEE", "#CC6677",  "#AA4499",   "#999933", "#882255", "#984EA3")
n = 10000
#relativePath = "C:/Users/lxy75/OneDrive/Documents/R-Project/mtDNA/Data/1130Prelim"
#relativePath_save = "C:/Users/lxy75/OneDrive/Documents/R-Project/mtDNA/Analysis/1130Prelim"

relativePath = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/20240213PowerRunFull/c17p6"
relativePath_save = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Analysis/20240320AlexPower"

# read all the meanDiffLL of the loops for this condition into a vector

allMeanDiffLL <- c()
for (i in 1:20){
      env <- new.env()
      load(paste0(relativePath,"/","loop",i,"/modelSmr.Rdata"), envir = env)
      allMeanDiffLL[i] <- env$smr2$Minus2LogLikelihood - env$smr1$Minus2LogLikelihood
}
library(psych)
describe(allMeanDiffLL)
plot(allMeanDiffLL)
PN <- ncol(env[["smr1"]][["dataSummary"]][["fam1"]])*length(env[["smr1"]][["dataSummary"]])
lamdaUnit <- mean(allMeanDiffLL)/PN
SSize <- 1: n
LamdaVec <- lamdaUnit*SSize
powVec <- as.numeric(lapply(LamdaVec,powerCal, df = 2))
df_power_unit <- data.frame(Nped = SSize, power = powVec)

