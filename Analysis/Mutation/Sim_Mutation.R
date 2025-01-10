# This script is designed for analyzing the influence of mtDNA mutations on the mtDNA effects estimate

# Load all the functions and csvs
source("/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/InitialDataPrep.R")


# we should use pedigree 12
ped = simulatePedigree(kpc = df_ped$k[12],
                       Ngen = df_ped$G[12],
                       sexR = df_ped$p[12],
                       marR = df_ped$r[12],
                       rd_kpc = TRUE)

# mutation designs:
# if the mutation happens, the mt2 matrix will work like seperate the descendants as its own pedigree.
# all other matrices will work just as normal. So we can just replace the mt2 matrix with the mutated mt2 matrix.
# Then we simulatee data based on the mutated matrices, but fit the model specified by the original matrices.
# conditions: happen at 2nd/3rd/4th generation
# conditions: happen 1/2/3 times

# # make a function to seperate the pedigree: drop a person from his/her parents
# #' @param ped a pedigree simulated from simulatePedigree function or the same format
# #' @param ID_drop the ID of the person to be dropped from his/her parents. 
# #' @param gen_drop the generation in which the randomly dropped person is. Will work if ID_drop is not specified.
# #' @param sex_drop the biological sex of the randomly dropped person.
# #' @param n_drop the number of times the mutation happens.
# dropLink = function(ped, 
#                     ID_drop = NA_integer_,
#                     gen_drop = 2,
#                     sex_drop = NA_character_,
#                     n_drop = 1){ 
#                         # a supporting function 
#                         resample = function(x, ...){x[sample.int(length(x), ...)]} 
#                          # check if the ID_drop is specified
#                         if(is.na(ID_drop)){
#                             # check if the sex_drop is specified
#                             if(is.na(sex_drop)){
#                                 ID_drop = resample(ped$ID[ped$gen == gen_drop & !is.na(ped$dadID) & !is.na(ped$momID)], n_drop)
#                             } else {
#                                 ID_drop = resample(ped$ID[ped$gen == gen_drop & !is.na(ped$dadID) & !is.na(ped$momID) & ped$sex == sex_drop], n_drop)
#                             }
#                             ped[ped$ID %in% ID_drop, c("dadID", "momID")] = NA_integer_
#                         } else {
#                             ped[ped$ID == ID_drop, c("dadID", "momID")] = NA_integer_
#                         }
#                         return(ped)
#                     }

# check if the function works. STATUS: GREEN. Function put in Functions
# ped_drop = dropLink(ped, gen_drop = 2, n_drop = 1)
# PlotPedigree(ped_drop)

# simulate and fit the model
# one mutation
# 2nd generation
# female missing
path_results = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/Mutation"

RunSim_Mutation = function(Var, kpc, Ngen, sexR, marR, n=10000, path_results, gen_drop, sex_drop, n_drop){
    library(OpenMx)
    library(mvtnorm)
    #Var Comb
    ad2 <- Var[[2]]
    cn2 <- Var[[3]]
    ce2 <- Var[[4]]
    mt2 <- Var[[5]]
    dd2 <- Var[[6]]
    am2 <- Var[[7]]
    ee2 <- Var[[8]]
    numfam <- round(n/famSizeCal(kpc,Ngen,marR))
    p.list <- list()
    p.list_missing <- list()
    modList <- list()
    modNames <- paste0("fam", 1:numfam)
    modList2 <- list()
    modNames2 <- paste0("fam", 1:numfam)
    modList3 <- list()
    modNames3 <- paste0("fam", 1:numfam)
    dat <- data.frame()
    for(i in 1: numfam){
        #set.seed(i*11)
        # Generate data and covariance matrices for the complete pedigree
        ped_temp <- simulatePedigree(kpc = kpc,
                                     Ngen = Ngen,
                                     sexR = sexR,
                                     marR = marR,
                                     rd_kpc = TRUE)
        #print(ped_temp)
        assign(paste0("ped",i),
               ped_temp)
        p.list[[i]] <- ped_temp
        names(p.list)[i] <- paste0("ped",i)
        
        Addmat <- as.matrix(ped2add(p.list[[i]], verbose = FALSE))
        Nucmat <- ped2cn(p.list[[i]])
        Extmat <- ped2ce(p.list[[i]])
        Mtdmat <- ped2mt_v3(p.list[[i]])
        Envmat <- diag(1,nrow = nrow(Addmat))
        dimnames(Envmat) <- dimnames(Addmat)
        Amimat <- Addmat*Mtdmat
        Dmgmat <- Addmat*Addmat
        #sumCov <- numeric()

        # Generate covariance matrices for the missing pedigree
        ped_drop = dropLink(ped_temp, gen_drop = gen_drop, sex_drop = sex_drop, n_drop = n_drop)
        assign(paste0("ped_drop",i),
               ped_drop)
        p.list_missing[[i]] <- ped_drop
        names(p.list_missing)[i] <- paste0("ped_drop",i)
        Addmat_missing <- as.matrix(ped2add(p.list[[i]], verbose = FALSE))
        Nucmat_missing <- ped2cn(p.list[[i]])
        Extmat_missing <- ped2ce(p.list[[i]])
        Mtdmat_missing <- ped2mt_v3(p.list_missing[[i]])
        Envmat_missing <- diag(1,nrow = nrow(Addmat_missing))
        dimnames(Envmat_missing) <- dimnames(Addmat_missing)
        Amimat_missing <- Addmat_missing*Mtdmat_missing
        Dmgmat_missing <- Addmat_missing*Addmat_missing
    
        ## generate data
        sumCov <- ad2*Addmat_missing + 
                  dd2*Addmat_missing*Addmat_missing + 
                  cn2*Nucmat_missing + 
                  ce2*Extmat_missing + 
                  mt2*Mtdmat_missing + 
                  am2*Addmat*Mtdmat_missing + 
                  ee2*Envmat_missing
        #set.seed(i*11)
        temp <- rmvnorm(1, sigma = sumCov)
        dat[i,1:ncol(Addmat)] <- as.double(temp)

        ytemp <- paste('S', rownames (Addmat))
        fsize <- nrow(Addmat)
        modList[[i]] <- mxModel(name=modNames[i],
                                   mxMatrix("Iden", nrow=fsize, ncol=fsize, name="I"), 
                                   mxMatrix("Unit", nrow=fsize, ncol=fsize, name='U'),
                                   mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Addmat, name="A"), 
                                   #mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Dmgmat, name="D"), 
                                   mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Nucmat, name="Cn"), 
                                   #mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Extmat, name="Ce"), 
                                   mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Amimat, name="Am"), 
                                   mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Mtdmat, name="Mt"),
                                   mxData(observed = matrix(temp, nrow=1, dimnames=list(NULL, ytemp)), type="raw", sort=FALSE),
                                   mxMatrix('Full', nrow=1, ncol=fsize, name='M', free=TRUE, labels='meanLI',
                                            dimnames=list(NULL, ytemp)),
                                   mxAlgebra ((A %x% ModelOne.Vad) 
                                              #+ (D %x% ModelOne.Vdd) 
                                              + (Cn %x% ModelOne.Vcn) 
                                              + (U %x% ModelOne.Vce) 
                                              + (Mt %x% ModelOne.Vmt) 
                                              + (Am %x% ModelOne.Vam) 
                                              + (I %x% ModelOne.Ver), 
                                              name="V", dimnames=list(ytemp, ytemp)),
                                   mxExpectationNormal(covariance='V', means='M'), 
                                   mxFitFunctionML()
        )
        
        
        ytemp2 <- paste('S', rownames (Addmat_missing))
        fsize2 <- nrow(Addmat_missing)
        modList2[[i]] <- mxModel(name=modNames2[i],
                                     mxMatrix("Iden", nrow=fsize2, ncol=fsize2, name="I"), 
                                     mxMatrix("Unit", nrow=fsize2, ncol=fsize2, name='U'),
                                     mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Addmat_missing, name="A"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Dmgmat, name="D"), 
                                     mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Nucmat_missing, name="Cn"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Extmat, name="Ce"), 
                                     mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Amimat_missing, name="Am"), 
                                     mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Mtdmat_missing, name="Mt"),
                                     mxData(observed = matrix(temp, nrow=1, dimnames=list(NULL, ytemp2)), type="raw", sort=FALSE),
                                     mxMatrix('Full', nrow=1, ncol=fsize2, name='M', free=TRUE, labels='meanLI',
                                              dimnames=list(NULL, ytemp2)),
                                     mxAlgebra ((A %x% ModelTwo.Vad) 
                                                #+ (D %x% ModelThree.Vdd) 
                                                + (Cn %x% ModelTwo.Vcn) 
                                                + (U %x% ModelTwo.Vce) 
                                                + (Mt %x% ModelTwo.Vmt) 
                                                + (Am %x% ModelTwo.Vam) 
                                                + (I %x% ModelTwo.Ver), 
                                                name="V", dimnames=list(ytemp2, ytemp2)),
                                     mxExpectationNormal(covariance='V', means='M'), 
                                     mxFitFunctionML()
        )

    }
    
    totalVar <- 1
    totalMea <- 0
    
    ObjectsKeep <- as.character(ls())
    # "l_ped","df_var","i","j","target_folder"
    ## the misspecified model
    Model1a <- mxModel(
        "ModelOne",
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ad2*totalVar, labels = "vad", name = "Vad", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = dd2*totalVar, labels = "vdd", name = "Vdd", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = cn2*totalVar, labels = "vcn", name = "Vcn", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ce2*totalVar, labels = "vce", name = "Vce", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = mt2*totalVar, labels = "vmt", name = "Vmt", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = am2*totalVar, labels = "vam", name = "Vam", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ee2*totalVar, labels = "ver", name = "Ver", lbound = 1e-10)
    )
    
    container <- mxModel('Model1b', Model1a, modList, mxFitFunctionMultigroup(modNames))
    container <- mxOption(container, 'Checkpoint Units', 'minutes')
    container <- mxOption(container, 'Checkpoint Count', 1)
    containerRun <- mxRun(container, intervals=FALSE, checkpoint=TRUE) 
    
    smr1 <- summary(containerRun)
    save(list = ls(envir = environment()), file = paste0(path_results,"/model1.RData"), envir = environment())
    rm(list = setdiff(ls(), c(ObjectsKeep, "ObjectsKeep", "smr1")))
    
    #### correct model
    
    Model2a <- mxModel(
        "ModelTwo",
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ad2*totalVar, labels = "vad", name = "Vad", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = dd2*totalVar, labels = "vdd", name = "Vdd", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = cn2*totalVar, labels = "vcn", name = "Vcn", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ce2*totalVar, labels = "vce", name = "Vce", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = mt2*totalVar, labels = "vmt", name = "Vmt", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = am2*totalVar, labels = "vam", name = "Vam", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ee2*totalVar, labels = "ver", name = "Ver", lbound = 1e-10)
    )
    
    container2 <- mxModel('Model2b', Model2a, modList2, mxFitFunctionMultigroup(modNames2))
    container2 <- mxOption(container2, 'Checkpoint Units', 'minutes')
    container2 <- mxOption(container2, 'Checkpoint Count', 1)
    containerRun2 <- mxRun(container2, intervals=FALSE, checkpoint=TRUE) 
    smr2 <- summary(containerRun2)
    save(list = ls(envir = environment()), file = paste0(path_results,"/model2.RData"), envir = environment())
    rm(list = setdiff(ls(), c(ObjectsKeep, "ObjectsKeep","path_results","smr1", "smr2")))
    
    #### save the summaries
    rm(list = setdiff(ls(), c("path_results","smr1", "smr2")))
     
     save(list = ls(envir = environment()), file = paste0(path_results,"/modelSmr.RData"), envir = environment())
 
     rm(list = ls())
}

RunSim_MissingLink = function(Var, kpc, Ngen, sexR, marR, n=10000, path_results, gen_drop, sex_drop, n_drop){
    library(OpenMx)
    library(mvtnorm)
    #Var Comb
    ad2 <- Var[[2]]
    cn2 <- Var[[3]]
    ce2 <- Var[[4]]
    mt2 <- Var[[5]]
    dd2 <- Var[[6]]
    am2 <- Var[[7]]
    ee2 <- Var[[8]]
    numfam <- round(n/famSizeCal(kpc,Ngen,marR))
    p.list <- list()
    p.list_missing <- list()
    modList <- list()
    modNames <- paste0("fam", 1:numfam)
    modList2 <- list()
    modNames2 <- paste0("fam", 1:numfam)
    modList3 <- list()
    modNames3 <- paste0("fam", 1:numfam)
    dat <- data.frame()
    for(i in 1: numfam){
        #set.seed(i*11)
        # Generate data and covariance matrices for the complete pedigree
        ped_temp <- simulatePedigree(kpc = kpc,
                                     Ngen = Ngen,
                                     sexR = sexR,
                                     marR = marR,
                                     rd_kpc = TRUE)
        #print(ped_temp)
        assign(paste0("ped",i),
               ped_temp)
        p.list[[i]] <- ped_temp
        names(p.list)[i] <- paste0("ped",i)
        
        Addmat <- as.matrix(ped2add(p.list[[i]], verbose = FALSE))
        Nucmat <- ped2cn(p.list[[i]])
        Extmat <- ped2ce(p.list[[i]])
        Mtdmat <- ped2mt_v3(p.list[[i]])
        Envmat <- diag(1,nrow = nrow(Addmat))
        dimnames(Envmat) <- dimnames(Addmat)
        Amimat <- Addmat*Mtdmat
        Dmgmat <- Addmat*Addmat
        #sumCov <- numeric()

        # Generate covariance matrices for the missing pedigree
        ped_drop = dropLink(ped_temp, gen_drop = gen_drop, sex_drop = sex_drop, n_drop = n_drop)
        assign(paste0("ped_drop",i),
               ped_drop)
        p.list_missing[[i]] <- ped_drop
        names(p.list_missing)[i] <- paste0("ped_drop",i)
        Addmat_missing <- as.matrix(ped2add(p.list_missing[[i]], verbose = FALSE))
        Nucmat_missing <- ped2cn(p.list_missing[[i]])
        Extmat_missing <- ped2ce(p.list_missing[[i]])
        Mtdmat_missing <- ped2mt_v3(p.list_missing[[i]])
        Envmat_missing <- diag(1,nrow = nrow(Addmat_missing))
        dimnames(Envmat_missing) <- dimnames(Addmat_missing)
        Amimat_missing <- Addmat_missing*Mtdmat_missing
        Dmgmat_missing <- Addmat_missing*Addmat_missing
    
        ## generate data
        sumCov <- ad2*Addmat_missing + 
                  dd2*Addmat_missing*Addmat_missing + 
                  cn2*Nucmat_missing + 
                  ce2*Extmat_missing + 
                  mt2*Mtdmat_missing + 
                  am2*Addmat*Mtdmat_missing + 
                  ee2*Envmat_missing
        #set.seed(i*11)
        temp <- rmvnorm(1, sigma = sumCov)
        dat[i,1:ncol(Addmat)] <- as.double(temp)

        ytemp <- paste('S', rownames (Addmat))
        fsize <- nrow(Addmat)
        modList[[i]] <- mxModel(name=modNames[i],
                                   mxMatrix("Iden", nrow=fsize, ncol=fsize, name="I"), 
                                   mxMatrix("Unit", nrow=fsize, ncol=fsize, name='U'),
                                   mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Addmat, name="A"), 
                                   #mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Dmgmat, name="D"), 
                                   mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Nucmat, name="Cn"), 
                                   #mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Extmat, name="Ce"), 
                                   mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Amimat, name="Am"), 
                                   mxMatrix("Symm", nrow=fsize, ncol=fsize, values=Mtdmat, name="Mt"),
                                   mxData(observed = matrix(temp, nrow=1, dimnames=list(NULL, ytemp)), type="raw", sort=FALSE),
                                   mxMatrix('Full', nrow=1, ncol=fsize, name='M', free=TRUE, labels='meanLI',
                                            dimnames=list(NULL, ytemp)),
                                   mxAlgebra ((A %x% ModelOne.Vad) 
                                              #+ (D %x% ModelOne.Vdd) 
                                              + (Cn %x% ModelOne.Vcn) 
                                              + (U %x% ModelOne.Vce) 
                                              + (Mt %x% ModelOne.Vmt) 
                                              + (Am %x% ModelOne.Vam) 
                                              + (I %x% ModelOne.Ver), 
                                              name="V", dimnames=list(ytemp, ytemp)),
                                   mxExpectationNormal(covariance='V', means='M'), 
                                   mxFitFunctionML()
        )
        
        
        ytemp2 <- paste('S', rownames (Addmat_missing))
        fsize2 <- nrow(Addmat_missing)
        modList2[[i]] <- mxModel(name=modNames2[i],
                                     mxMatrix("Iden", nrow=fsize2, ncol=fsize2, name="I"), 
                                     mxMatrix("Unit", nrow=fsize2, ncol=fsize2, name='U'),
                                     mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Addmat_missing, name="A"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Dmgmat, name="D"), 
                                     mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Nucmat_missing, name="Cn"), 
                                     #mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Extmat, name="Ce"), 
                                     mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Amimat_missing, name="Am"), 
                                     mxMatrix("Symm", nrow=fsize2, ncol=fsize2, values=Mtdmat_missing, name="Mt"),
                                     mxData(observed = matrix(temp, nrow=1, dimnames=list(NULL, ytemp2)), type="raw", sort=FALSE),
                                     mxMatrix('Full', nrow=1, ncol=fsize2, name='M', free=TRUE, labels='meanLI',
                                              dimnames=list(NULL, ytemp2)),
                                     mxAlgebra ((A %x% ModelTwo.Vad) 
                                                #+ (D %x% ModelThree.Vdd) 
                                                + (Cn %x% ModelTwo.Vcn) 
                                                + (U %x% ModelTwo.Vce) 
                                                + (Mt %x% ModelTwo.Vmt) 
                                                + (Am %x% ModelTwo.Vam) 
                                                + (I %x% ModelTwo.Ver), 
                                                name="V", dimnames=list(ytemp2, ytemp2)),
                                     mxExpectationNormal(covariance='V', means='M'), 
                                     mxFitFunctionML()
        )

    }
    
    totalVar <- 1
    totalMea <- 0
    
    ObjectsKeep <- as.character(ls())
    # "l_ped","df_var","i","j","target_folder"
    ## the misspecified model
    Model1a <- mxModel(
        "ModelOne",
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ad2*totalVar, labels = "vad", name = "Vad", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = dd2*totalVar, labels = "vdd", name = "Vdd", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = cn2*totalVar, labels = "vcn", name = "Vcn", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ce2*totalVar, labels = "vce", name = "Vce", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = mt2*totalVar, labels = "vmt", name = "Vmt", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = am2*totalVar, labels = "vam", name = "Vam", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ee2*totalVar, labels = "ver", name = "Ver", lbound = 1e-10)
    )
    
    container <- mxModel('Model1b', Model1a, modList, mxFitFunctionMultigroup(modNames))
    container <- mxOption(container, 'Checkpoint Units', 'minutes')
    container <- mxOption(container, 'Checkpoint Count', 1)
    containerRun <- mxRun(container, intervals=FALSE, checkpoint=TRUE) 
    
    smr1 <- summary(containerRun)
    save(list = ls(envir = environment()), file = paste0(path_results,"/model1.RData"), envir = environment())
    rm(list = setdiff(ls(), c(ObjectsKeep, "ObjectsKeep", "smr1")))
    
    #### correct model
    
    Model2a <- mxModel(
        "ModelTwo",
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ad2*totalVar, labels = "vad", name = "Vad", lbound = 1e-10),
        #mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = dd2*totalVar, labels = "vdd", name = "Vdd", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = cn2*totalVar, labels = "vcn", name = "Vcn", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ce2*totalVar, labels = "vce", name = "Vce", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = mt2*totalVar, labels = "vmt", name = "Vmt", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = am2*totalVar, labels = "vam", name = "Vam", lbound = 1e-10),
        mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = ee2*totalVar, labels = "ver", name = "Ver", lbound = 1e-10)
    )
    
    container2 <- mxModel('Model2b', Model2a, modList2, mxFitFunctionMultigroup(modNames2))
    container2 <- mxOption(container2, 'Checkpoint Units', 'minutes')
    container2 <- mxOption(container2, 'Checkpoint Count', 1)
    containerRun2 <- mxRun(container2, intervals=FALSE, checkpoint=TRUE) 
    smr2 <- summary(containerRun2)
    save(list = ls(envir = environment()), file = paste0(path_results,"/model2.RData"), envir = environment())
    rm(list = setdiff(ls(), c(ObjectsKeep, "ObjectsKeep","path_results","smr1", "smr2")))
    
    #### save the summaries
    rm(list = setdiff(ls(), c("path_results","smr1", "smr2")))
     
     save(list = ls(envir = environment()), file = paste0(path_results,"/modelSmr.RData"), envir = environment())
 
     rm(list = ls())
}