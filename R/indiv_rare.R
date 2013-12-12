#' Individual base rarefaction
#' @description Produce individual based rarefaction curves based on Colwell et al. 2012 (doi: 10.1093/jpe/rtr044),  allows for both interpolation and extrapolation
#' @param smat a vector of species abundances or a site by species matrix 
#' @param extrap the number of individuals to extrapolate out to
#' @details If you extrapolate too far the algorithm results may not make much sense.  Currently this method uses just the basic Chao1 estimatory, but eventually will include the ACE estimator. Just one value is accepted here, it will extend all samples in a matrix out to this common end point for comparison sake.
#' @return a vector of interpolated and extrapolated values if a vector is input, otherwise a list with sites arranged in the same order as they were in the input matrix
#' @examples \dontrun{
#' ## Data from Colwell 2012 table 4
#' data(costaRicaTrees)  
#' LEPrare <- indiv_rare(old_growth)
#' plot(LEPrare, type='l', xlab="Individuals",ylab="Species") 
#' 
#' ## Matrix example with data from Work et al 2010
#' data(work)
#' work_rare <- indiv_rare(as.matrix(work))
#' ### Plot results
#' plot(x = unlist(lapply(work_rare,length)),unlist(lapply(work_rare,max)),ylim=c(0,max(unlist(work_rare))),xlim=c(0,max(unlist(lapply(work_rare,length)))), ylab = "Species", xlab ="Individuals")
#' for(i in 1:length(work_rare)){
#' lines(work_rare[[i]])
#' }
#' 
#' }
#' @export indiv_rare



indiv_rare <- function(smat,extrap = NULL){
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




