\name{indiv_rare}
\alias{indiv_rare}
\title{Individual base rarefaction}
\usage{
  indiv_rare(smat)
}
\arguments{
  \item{smat}{a vector of species abundances or a site by
  species matrix}
}
\value{
  a vector of interpolated and extrapolated values if a
  vector is input, otherwise a list with sites arranged in
  the same order as they were in the input matrix
}
\description{
  Produce individual based rarefaction curves based on
  Colwell et al. 2012 (doi: 10.1093/jpe/rtr044), allows for
  both interpolation and extrapolation
}
\examples{
\dontrun{
## Data from Colwell 2012 table 4
ni <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 16, 18, 19, 20, 25, 38, 39, 40, 46, 52, 55)
count <- c(46, 30, 16, 12, 6, 5, 3, 4, 5, 4, 1, 3, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1)
LEP_cost_og <- vector()
for(i in 1:length(ni)){
 LEP_cost_og <- c(LEP_cost_og,rep(ni[i],count[i]))}

LEPrare <- indiv_rare(LEP_cost_og)
plot(LEPrare, type='l', xlab="Individuals",ylab="Species")

## Matrix example
emend <- as.data.frame(read.csv("http://www.jennajacobs.org/R/EMEND.csv", row.names=1))
emend_rare <- indiv_rare(as.matrix(emend))
### Plot results
plot(x = unlist(lapply(emend_rare,length)),unlist(lapply(emend_rare,max)),ylim=c(0,max(unlist(emend_rare))),xlim=c(0,max(unlist(lapply(emend_rare,length)))), ylab = "Species", xlab ="Individuals")
for(i in 1:length(emend_rare)){
lines(emend_rare[[i]])
}

}
}

