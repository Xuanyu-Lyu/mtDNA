# a power summary for all different conditions
relativePath = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/20240213PowerRunFull"
relativePath_save = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Analysis/20240321PowerMean50"
powerCal <- function(lamda,df){
      1- pchisq(qchisq(1-.05, df), df, lamda)
      
}
n = 10000

getPowerDf <- function(condition, n = 20000){
        # read all the meanDiffLL of the loops for this condition into a vector
    powerCal <- function(lamda,df){
      1- pchisq(qchisq(1-.05, df), df, lamda)
      
    }
    allMeanDiffLL <- c()
    for (i in 1:100){
        # read all the meanDiffLL of the loops for this condition into a vector

    
        env <- new.env()
        load(paste0(relativePath,condition,"loop",i,"/modelSmr.Rdata"), envir = env)
        allMeanDiffLL[i] <- env$smr2$Minus2LogLikelihood - env$smr1$Minus2LogLikelihood
    }

    PN <- ncol(env[["smr1"]][["dataSummary"]][["fam1"]])*length(env[["smr1"]][["dataSummary"]])
    meanDiffLL_mtam <- mean(allMeanDiffLL)
    lamdaUnit <- meanDiffLL_mtam/PN
    SSize <- 1: n
    LamdaVec <- lamdaUnit*SSize
    powVec <- as.numeric(lapply(LamdaVec,powerCal, df = 2))

    df <- data.frame(Nped = SSize, 
                        power = powVec, 
                        Combination = rep("power", n))
    return(df)
}



printSS <- function(df){
    # find the lowest Nped gives power > .4, .6, .8 and =1 respectively, and put them in a vector
    Nped_4 <- min(df[df$power > .4,]$Nped)
    Nped_6 <- min(df[df$power > .6,]$Nped)
    Nped_8 <- min(df[df$power > .8,]$Nped)
    Nped_1 <- min(df[df$power > .9999,]$Nped)
    return(c(Nped_4, Nped_6, Nped_8, Nped_1))
}

df_c1p4 <- getPowerDf("/c1p4/")
printSS(df_c1p4) 

df_c2p4 <- getPowerDf("/c2p4/")
printSS(df_c2p4) 

df_c3p4 <- getPowerDf("/c3p4/", n = 2e6)
printSS(df_c3p4) 

df_c4p4 <- getPowerDf("/c4p4/", n = 2e7)
printSS(df_c4p4) 



df_c1p1 <- getPowerDf("/c1p1/")
printSS(df_c1p1) 

df_c2p1 <- getPowerDf("/c2p1/")
printSS(df_c2p1) 

df_c3p1 <- getPowerDf("/c3p1/", n = 2e6)
printSS(df_c3p1) 

df_c4p1 <- getPowerDf("/c4p1/", n = 2e6)
printSS(df_c4p1)


df_c2p1 <- getPowerDf("/c2p1/")
printSS(df_c2p1) 
df_c2p4 <- getPowerDf("/c2p4/")
printSS(df_c2p4)
df_c2p5 <- getPowerDf("/c2p5/")
printSS(df_c2p5)

df_c2p1 <- getPowerDf("/c2p1/")
printSS(df_c2p1) 
df_c2p2 <- getPowerDf("/c2p2/", n = 30000)
printSS(df_c2p2) 
df_c2p3 <- getPowerDf("/c2p3/", n = 30000)
printSS(df_c2p3) 

(19396-15182)/15182

install.packages("BGmisc")
library(BGmisc)
# Create a new PNG file
png(filename = "pedigree_plot1.png", width = 800, height = 800, res = 200)

ped1 <- simulatePedigree(Ngen = 3, kpc = 4, marR = .7,rd_kpc = TRUE)

# Generate the pedigree plot
x1 = plotPedigree(ped1)

# Close the PNG file
dev.off()
# save the pedigree plot



