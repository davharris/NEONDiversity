#' Calculate Hurlbert's PIE 
#' @description Calculate Hurlbert's PIE (Probability of an Interspecific Encounter)
#' @param sp_list a vector of species abundances, or a site (rows) x species (columns) matrix
#' @details This will calculate Hurlbert's PIE after Hurlbert 1971.  The equation used is \eqn{PIE = \begin{bmatrix} \frac{N}{N-1} \end{bmatrix} \begin{bmatrix} 1 - \sum_{i}^{N}\begin{pmatrix} \frac{N_i}{N} \end{pmatrix}^2 \end{bmatrix}}{ PIE = \begin{bmatrix} \frac{N}{N-1} \end{bmatrix} \begin{bmatrix} 1 - \sum_{i}^{N}\begin{pmatrix} \frac{N_i}{N} \end{pmatrix}^2 \end{bmatrix} }
#' @return a scalar (for a vector of input) or vector (for matrix input) of PIE values
#' @examples \dontrun{
#' ## Generate some data
#' fake_dat <- matrix(rpois(100,100),ncol=10,nrow=10)
#' 
#' ### Return 10 values for PIE
#' 
#' h_pie(fake_dat)
#' 
#' }
#' 
#' @export h_pie


h_pie <- function(sp_list){
  
  #convert vector to matrix
  if(is.vector(sp_list)){
    sp_list <- t(as.matrix(sp_list))
  }
  
  out <- apply(sp_list,1, NEONDiversity:::pie_calc)
  return(out)
  
}

pie_calc <- function(vec) {
  
  N <- sum(vec)
  a <- (N / ( N -1 ))
  b <- 1 - sum((vec/N)^2)
  return(a*b)
  
}

