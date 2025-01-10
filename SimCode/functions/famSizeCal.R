# A function to return a vector for family sizes in each generation given:
# Parameters:
# kpc: kids per couple 
# Ngen: number of generations
# marR: marriage rate
allGens <- function(kpc, Ngen, marR){
   if (Ngen < 1){
      stop("The number of generations should be integers greater or equal than 1")
   }
   if (Ngen == 1){
      allGens = 2
   }
   if (Ngen >= 2){
      allGens = sizeAllGens(kpc=kpc, Ngen = Ngen, marR = marR)
   } else{
      stop()
   }
   return(allGens)
}


# A function to calculate the expected family size given 
# 1) Kids per couple
# 2) Number of generations
# 3) Marriage rate


# Parameters:
# kpc: kids per couple 
# Ngen: number of generations
# marR: marriage rate
famSizeCal <- function(kpc, Ngen, marR){
   if (Ngen < 1){
      stop("The number of generations should be integers greater or equal than 1")
   }
   if (Ngen == 1){
      size = 2
   }
   if (Ngen >= 2){
      allGens = sizeAllGens(kpc=kpc, Ngen = Ngen, marR = marR)
      size = sum(allGens)
   } else{
      stop()
   }
   return(size)
}

# A function to get a vector of of family members in all generations for all  using an exponential function of kpc, Ngen and marR.
# Ideally, all there parameters should be inherited from the parent function famSizeCal
sizeAllGens <- function(kpc, Ngen, marR){
   Nmid <- Ngen -2
   midGens <- numeric(length = Nmid)
   
   for(i in 2 : (Ngen - 1)){
      midGens[i-1] <- kpc^(i-1) * marR^(i-2) * (1+marR)
      midGens[i-1] <- ceiling(midGens[i-1])
   }
   
   lastGen <- ceiling(kpc^(Ngen-1) * marR^(Ngen-2))
   allGens <- c(2,midGens,lastGen)
   #print(allGens)
   return(allGens)
}

# test the function: sizeMidCal green
#results <- sizeAllGens(5,6,1)
#results <- allGens(3,4,2/3)
#results <- famSizeCal(3,4,2/3)
# test the function famSizeCal: g