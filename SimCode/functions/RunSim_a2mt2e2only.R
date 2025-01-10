# A function to run a specific simulation condition
# Var: A vector; selected variance component
# Ped: A dataframe; selected pedigree structure
# data_ini: a file path-the initializing Rdata file containing the pedigrees and variance combinations for simulation
# path_results: a path of the folder indicating the place for storing the simulation results.
RunSim_a2mt2e2only <- function(Var, Ped,n=10000, path_results){
    source("~/R-Project/mtDNA_mt2/Functions/HelperFunctions.R")
    library(OpenMx)
    library(mvtnorm)
    
    Addmat <- as.matrix(ped2add(Ped, verbose = TRUE))
    Nucmat <- ped2cn(Ped)
    Extmat <- ped2ce(Ped)
    Mtdmat <- ped2mt_v3(Ped)
    Envmat <- diag(1,nrow = nrow(Addmat))
    dimnames(Envmat) <- dimnames(Addmat)
    Amimat <- Addmat*Mtdmat
    Dmgmat <- Addmat*Addmat
    
    #Var Comb
    ad2 <- Var[[2]]
    #print(class(ad2))
    cn2 <- Var[[3]]
    ce2 <- Var[[4]]
    mt2 <- Var[[5]]
    dd2 <- Var[[6]]
    am2 <- Var[[7]]
    ee2 <- Var[[8]]
    
    ## generate data
    sumCov <- ad2*Addmat + dd2*Addmat*Addmat + cn2*Nucmat + ce2*Extmat + mt2*Mtdmat + am2*Addmat*Mtdmat + ee2*Envmat
    set.seed(142)
    numfam <- round(n/nrow(Addmat))
    #print(class(numfam))
    dat <- rmvnorm(numfam, sigma = sumCov)
    
    totalVar <- 1
    totalMea <- 0
    
    ObjectsKeep <- as.character(ls())
    # "l_ped","df_var","i","j","target_folder"
    ## the full model
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
    
    save(list = ls(envir = environment()), file = paste0(path_results,"/model1.RData"), envir = environment())
    
    #save.image(file = paste0(path_results,"/model1.RData"))
    
    rm(list = setdiff(ls(), c(ObjectsKeep, "ObjectsKeep", "smr1")))
    
    ##### Model excluding mt and am
    
    Model2a <- mxModel(
        "ModelThree",
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ad2*totalVar, labels = "vad", name = "Vad", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = dd2*totalVar, labels = "vdd", name = "Vdd", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = cn2*totalVar, labels = "vcn", name = "Vcn", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ce2*totalVar, labels = "vce", name = "Vce", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = mt2*totalVar, labels = "vmt", name = "Vmt", lbound = 1e-10),
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
                                     mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Addmat, name="A"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Dmgmat, name="D"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Nucmat, name="Cn"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Extmat, name="Ce"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Amimat, name="Am"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Mtdmat, name="Mt"),
                                     mxData(observed = matrix(ll2[[afam2]], nrow=1, dimnames=list(NULL, ytemp2)), type="raw", sort=FALSE),
                                     mxMatrix('Full', nrow=1, ncol=fsize2, name='M', free=TRUE, labels='meanLI',
                                              dimnames=list(NULL, ytemp2)),
                                     mxAlgebra ((A %x% ModelThree.Vad) 
                                                #+ (D %x% ModelThree.Vdd) 
                                               # + (Cn %x% ModelThree.Vcn) 
                                                #+ (U %x% ModelThree.Vce) 
                                                #+ (Mt %x% ModelThree.Vmt) 
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
    save(list = ls(envir = environment()), file = paste0(path_results,"/model2.RData"), envir = environment())
    
    #rm(list = setdiff(ls(), c(ObjectsKeep,"ObjectsKeep", "smr1", "smr2")))
    
    rm(list = setdiff(ls(), c("path_results","smr1", "smr2")))
    
    save(list = ls(envir = environment()), file = paste0(path_results,"/modelSmr.RData"), envir = environment())
    #save.image(file = paste0(path_results,"/modelSmr.RData") )
    rm(list = ls())
}