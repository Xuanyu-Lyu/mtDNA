#devtools::install_github("R-Computing-Lab/BGmisc")
## Load the package
library(BGmisc)
library(OpenMx)
library(mvtnorm)

## k 1-100 Ped1 = simulatePedigree(kpc = 4, Ngen = 5, sexR = .5, marR = 1)


for (k in 1:100) {
    #Set the wd
    path_results <- "~/R-Project/mtDNA_mt2/MissLink"
    
    # Simulate the pedigree first
    set.seed(63)
    Ped1 = simulatePedigree(kpc = 4, Ngen = 5, sexR = .5, marR = 1)
    
    plotPedigree(Ped1)
    
    # Drop a kin link in the pedigree
    Ped1miss <- Ped1
    Ped1miss[Ped1miss["ID"]==10024, c("dadID", "momID")] <- c(NA,NA)
    
    plotPedigree(Ped1miss)
    # Generate related matrices for the complete pedigree
    Addmat <- as.matrix(ped2add(Ped1, verbose = TRUE))
    Nucmat <- ped2cn(Ped1)
    Extmat <- ped2ce(Ped1)
    Mtdmat <- ped2mit(Ped1)
    Envmat <- diag(1,nrow = nrow(Addmat))
    dimnames(Envmat) <- dimnames(Addmat)
    Amimat <- Addmat*Mtdmat
    Dmgmat <- Addmat*Addmat
    
    # Generate related matrices for the complete pedigree
    Addmat_miss <- as.matrix(ped2add(Ped1miss, verbose = TRUE))
    Nucmat_miss <- ped2cn(Ped1miss)
    Extmat_miss <- ped2ce(Ped1miss)
    Mtdmat_miss <- ped2mit(Ped1miss)
    Envmat_miss <- diag(1,nrow = nrow(Addmat_miss))
    dimnames(Envmat_miss) <- dimnames(Addmat_miss)
    Amimat_miss <- Addmat_miss*Mtdmat_miss
    Dmgmat_miss <- Addmat_miss*Addmat_miss
    
    
    #Var Comb
    df_var <- read.csv("~/R-Project/mtDNA_mt2/InitialData/VarianceComb_a2mt2.csv")
    rowR <- 1
    ad2 <- df_var[rowR,2]
    cn2 <- df_var[rowR,3]
    ce2 <- df_var[rowR,4]
    mt2 <- df_var[rowR,5]
    dd2 <- df_var[rowR,6]
    am2 <- df_var[rowR,7]
    ee2 <- df_var[rowR,8]
    sumCov <- ad2*Addmat + dd2*Addmat*Addmat + cn2*Nucmat + ce2*Extmat + mt2*Mtdmat + am2*Addmat*Mtdmat + ee2*Envmat
    
    n <- 5000
    numfam <- round(n/nrow(Addmat))
    set.seed(60+k)
    dat <- rmvnorm(numfam, sigma = sumCov)
    
    totalVar <- 1
    totalMea <- 0
    ObjectsKeep <- as.character(ls())
    
    Model1a <- mxModel(
        "ModelOne",
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ad2*totalVar, labels = "vad", name = "Vad", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = dd2*totalVar, labels = "vdd", name = "Vdd", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = cn2*totalVar, labels = "vcn", name = "Vcn", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ce2*totalVar, labels = "vce", name = "Vce", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = mt2*totalVar, labels = "vmt", name = "Vmt", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = am2*totalVar, labels = "vam", name = "Vam", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ee2*totalVar, labels = "ver", name = "Ver", lbound = 1e-10)
    )
    
    ll <- list()
    for (i in 1:numfam){
        ll[[i]] <- dat[i,]
    }
    
    modList <- list()
    modNames <- paste0("fam", 1:numfam)
    
    for(afam in 1:numfam){
        ytemp <- paste('S', rownames (Addmat))
        fsize <- nrow(Addmat)
        modList[[afam]] <- mxModel(name=modNames[afam],
                                   mxMatrix("Iden", nrow=fsize, ncol=fsize, name="I"), 
                                   #mxMatrix("Unit", nrow=fsize, ncol=fsize, name='U'),
                                   mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Addmat, name="A"), 
                                   #mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Dmgmat, name="D"), 
                                   #mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Nucmat, name="Cn"), 
                                   #mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Extmat, name="Ce"), 
                                   #mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Amimat, name="Am"), 
                                   mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Mtdmat, name="Mt"),
                                   mxData(observed = matrix(ll[[afam]], nrow=1, dimnames=list(NULL, ytemp)), type="raw", sort=FALSE),
                                   mxMatrix('Full', nrow=1, ncol=fsize, name='M', free=TRUE, labels='meanLI',
                                            dimnames=list(NULL, ytemp)),
                                   mxAlgebra ((A %x% ModelOne.Vad) 
                                              #+ (D %x% ModelOne.Vdd) 
                                              #+ (Cn %x% ModelOne.Vcn) 
                                              #+ (U %x% ModelOne.Vce) 
                                              + (Mt %x% ModelOne.Vmt) 
                                              #+ (Am %x% ModelOne.Vam) 
                                              + (I %x% ModelOne.Ver), 
                                              name="V", dimnames=list(ytemp, ytemp)),
                                   mxExpectationNormal(covariance='V', means='M'), 
                                   mxFitFunctionML()
        )
    }
    container <- mxModel('Model1b', Model1a, modList, mxFitFunctionMultigroup(modNames))
    container <- mxOption(container, 'Checkpoint Units', 'minutes')
    container <- mxOption(container, 'Checkpoint Count', 1)
    containerRun <- mxRun(container, intervals=FALSE, checkpoint=TRUE) 
    
    smr1 <- summary(containerRun)
    
    save(list = ls(envir = environment()), file = paste0(path_results,paste0("/model1_",k,".RData")), envir = environment())
    
    
    rm(list = setdiff(ls(), c(ObjectsKeep, "ObjectsKeep", "smr1")))
    
    ##### Model-fit the missing-link model with the correct data
    
    Model2a <- mxModel(
        "ModelThree",
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ad2*totalVar, labels = "vad", name = "Vad", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = dd2*totalVar, labels = "vdd", name = "Vdd", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = cn2*totalVar, labels = "vcn", name = "Vcn", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ce2*totalVar, labels = "vce", name = "Vce", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = mt2*totalVar, labels = "vmt", name = "Vmt", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = am2*totalVar, labels = "vam", name = "Vam", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ee2*totalVar, labels = "ver", name = "Ver", lbound = 1e-10)
    )
    
    ll2 <- list()
    for (i in 1:numfam){
        ll2[[i]] <- dat[i,]
    }
    
    modList2 <- list()
    modNames2 <- paste0("fam", 1:numfam)
    
    for(afam2 in 1:numfam){
        ytemp2 <- paste('S', rownames (Addmat))
        fsize2 <- nrow(Addmat)
        modList2[[afam2]] <- mxModel(name=modNames2[afam2],
                                     mxMatrix("Iden", nrow=fsize2, ncol=fsize2, name="I"), 
                                     #mxMatrix("Unit", nrow=fsize2, ncol=fsize2, name='U'),
                                     mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Addmat_miss, name="A"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Dmgmat_miss, name="D"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Nucmat_miss, name="Cn"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Extmat_miss, name="Ce"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Amimat_miss, name="Am"), 
                                     mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Mtdmat_miss, name="Mt"),
                                     mxData(observed = matrix(ll2[[afam2]], nrow=1, dimnames=list(NULL, ytemp2)), type="raw", sort=FALSE),
                                     mxMatrix('Full', nrow=1, ncol=fsize2, name='M', free=TRUE, labels='meanLI',
                                              dimnames=list(NULL, ytemp2)),
                                     mxAlgebra ((A %x% ModelThree.Vad) 
                                                #+ (D %x% ModelThree.Vdd) 
                                                #+ (Cn %x% ModelThree.Vcn) 
                                                #+ (U %x% ModelThree.Vce) 
                                                + (Mt %x% ModelThree.Vmt) 
                                                #+ (Am %x% ModelThree.Vam) 
                                                + (I %x% ModelThree.Ver), 
                                                name="V", dimnames=list(ytemp2, ytemp2)),
                                     mxExpectationNormal(covariance='V', means='M'), 
                                     mxFitFunctionML()
        )
    }
    container2 <- mxModel('Model2b', Model2a, modList2, mxFitFunctionMultigroup(modNames2))
    container2 <- mxOption(container2, 'Checkpoint Units', 'minutes')
    container2 <- mxOption(container2, 'Checkpoint Count', 1)
    containerRun2 <- mxRun(container2, intervals=FALSE, checkpoint=TRUE) 
    smr2 <- summary(containerRun2)
    
    #save.image(file = paste0(path_results,"/model2.RData"))
    save(list = ls(envir = environment()), file = paste0(path_results,paste0("/model1_",k,".RData")), envir = environment())
    
    rm(list = setdiff(ls(), c(ObjectsKeep,"ObjectsKeep", "smr1", "smr2")))
    
    #rm(list = setdiff(ls(), c("path_results","smr1", "smr2")))
    
    save(list = ls(envir = environment()), file = paste0(path_results,paste0("/modelSmr_",k,".RData")), envir = environment())
    #save.image(file = paste0(path_results,"/modelSmr.RData") )
    rm(list = ls())
    
    # load("~/R-Project/mtDNA_mt2/MissLink/modelSmr.RData")
    # smr1$parameters
    # smr2$parameters
    
}
 