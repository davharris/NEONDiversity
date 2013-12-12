#' Chao1 estimator
#' @description Calculates the Chao1 estimator based on Colwell 2012 and Chao 1984
#' @param spp a vector or a site x species matrix
#' @param f1 a scalar or vector of the frequency of singletons see details
#' @param f2 a scalar or vector of the frequency of doubletons see details
#' @details The function can take either a site by species matrix or alternatively the parameters f1, f2 and S.  Those are calculated from the matrix / vector, but don't need to be if that is all the information you have
#' @return a vector of estimates unobserved species which can be added to S for an estimate of the total number of species
#' @import plyr
#' @examples \dontrun{
#' data(costaRicaTrees)
#' chao1(new_growth)
#' 
#' }
#' 
#' @export chao1

chao1 <- function(spp = NULL, f1 = NULL ,f2 = NULL){
  ## Type conversion
  if(is.data.frame(spp)){
    spp <- as.matrix(spp)
  }
  if(!is.null(spp) && is.vector(spp)){
    # get counts
    f1 <- sum(spp == 1)
    f2 <- sum(spp == 2)
  } else if(!is.null(spp) && is.matrix(spp)){
    f1 <- apply(spp,1,function(x) sum(x == 1))
    f2 <- apply(spp,1,function(x) sum(x == 2))
  }
  f0 <- rep(NA,length(f1))
 for(i in 1:length(f1)){
   if(f2[i] > 0){
     f0[i] <- f1[i]^2 / (2*f2[i])
   } else if(f2[i] == 0){
     f0[i] <- f1[i]*(f1[i] - 1) / (2*(f2[i] + 1))
   }
 }
  return(f0)
  
}