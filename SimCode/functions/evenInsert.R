# this is a function to insert m elements in a length n vector evenly. This function supports the balanced sex ratio between spouses
# @m: a vector 
# @n: a vector
evenInsert <- function(m,n){
      if (length(m) > length(n)){
            temp <- m
            m <- n
            n <- temp
      }
      
      #idx <- numeric()
      for (i in 1:length(m)){
            names(m)[i] <- ceiling(i*length(n)/length(m))
      }
      #print(m)
      
      names(n) <- 1:length(n)
      #print(n)
      
      vec <- c(m,n)
      vec <- vec[order(as.numeric(names(vec)))]
      vec <- unname(vec)
      
      return(vec)
}

# evenInsert(rep("a", 100),
#            rep("b",5))
