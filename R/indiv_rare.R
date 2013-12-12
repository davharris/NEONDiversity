#' Individual base rarefaction
#' @description Produce individual based rarefaction curves based on Colwell et al. 2012 (doi: 10.1093/jpe/rtr044),  allows for both interpolation and extrapolation
#' @param smat a vector of species abundances or a site by species matrix 
#' @return a vector of interpolated and extrapolated values if a vector is input, otherwise a list with sites arranged in the same order as they were in the input matrix
#' @examples \dontrun{
#' ## Data from Colwell 2012 table 4
#' ni <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 16, 18, 19, 20, 25, 38, 39, 40, 46, 52, 55)
#' count <- c(46, 30, 16, 12, 6, 5, 3, 4, 5, 4, 1, 3, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1)
#' LEP_cost_og <- vector()
#' for(i in 1:length(ni)){
#'  LEP_cost_og <- c(LEP_cost_og,rep(ni[i],count[i]))}
#'  
#' LEPrare <- indiv_rare(LEP_cost_og)
#' plot(LEPrare, type='l', xlab="Individuals",ylab="Species") 
#' 
#' ## Matrix example
#' emend <- as.data.frame(read.csv("http://www.jennajacobs.org/R/EMEND.csv", row.names=1))
#' emend_rare <- indiv_rare(as.matrix(emend))
#' ### Plot results
#' plot(x = unlist(lapply(emend_rare,length)),unlist(lapply(emend_rare,max)),ylim=c(0,max(unlist(emend_rare))),xlim=c(0,max(unlist(lapply(emend_rare,length)))), ylab = "Species", xlab ="Individuals")
#' for(i in 1:length(emend_rare)){
#' lines(emend_rare[[i]])
#' }
#' 
#' }
#' @export indiv_rare



indiv_rare <- function(smat){
  if(is.vector(smat)){
    return(NEONDiversity:::vec_rare(smat))
  }
  
  if(is.data.frame(smat)){
    smat <- as.matrix(smat)
  }
  
  if(is.matrix(smat)){
    out <- list()
    for(i in 1:dim(smat)[1]){
      out[[i]] <- NEONDiversity:::vec_rare(as.vector(smat[i,]))
    }
      return(out)
    
  }
  
}

#' Individual base rarefaction
#' @description calculated individual based rarefaction for a single site
#' @param spp a vector of species abundances
#' @return interpolated values for all possible individual subsamples

vec_rare <- function(spp) {

## Get variable names from Colwell et al 2012 eq 4  
sobs <- length(spp)
n <- sum(spp)
freqk <- table(spp)
kclass <- as.numeric(names(freqk))


m_out <- vector()
for(m in 1:n){
  tval <- sobs - NEONDiversity:::rhs(spp,n,m,kclass,freqk)
  m_out <- c(m_out,tval)
}
return(m_out)

}


#' Calculate alpha function
#' @description calculate the alpha function using log factorials from Colwell 2012 see eqn 4
#' @param n the total number of individuals in the sample
#' @param m the number of individuals to interpolate S for.
#' @param k the number of individuals in a given frequency class 
#' @return Alpha summed over all k 
alpha <- function(n,k,m) {
  if(k <= (n-m)){
    out <- exp((lfactorial((n-k)) + lfactorial((n-m))) - (lfactorial(n) + lfactorial((n-k-m))))
    return(out)
  }
  else{
    return(0)
  }
}

#' the right hand side of eqn 4 in Colwell 2012
#' @description calculate the right hand side of eqn 4, works by only calculating values where frequency classes exist and summing them
#' @param vec the species abundance vector
#' @param n the total number of individuals in the sample
#' @param m the number of individuals to interpolate S for.
#' @param kclass the number of individuals in a given frequency class 
#' @param freqk a vector of frequency classes 
#' @details all these values can be calculated from the vector of species in the initial function call but I pass them for transparency's sake
#' @return the full rhs of eqn 4 as a scalar
rhs <- function(vec,n,m, kclass,freqk) {
  out <- vector()
  for(k in 1:length(kclass)){

    out <- c(out,NEONDiversity:::alpha(n,kclass[k],m) * freqk[k]) 
  }
return(sum(out,na.rm=T))
}

