\name{indiv_rare}
\alias{indiv_rare}
\title{Individual base rarefaction}
\usage{
  indiv_rare(smat, extrap = NULL)
}
\arguments{
  \item{smat}{a vector of species abundances or a site by
  species matrix}

  \item{extrap}{the number of individuals to extrapolate
  out to}
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
\details{
  If you extrapolate too far the algorithm results may not
  make much sense.  Currently this method uses just the
  basic Chao1 estimatory, but eventually will include the
  ACE estimator. Just one value is accepted here, it will
  extend all samples in a matrix out to this common end
  point for comparison sake.
}
\examples{
\dontrun{
## Data from Colwell 2012 table 4
data(costaRicaTrees)
LEPrare <- indiv_rare(old_growth)
plot(LEPrare, type='l', xlab="Individuals",ylab="Species")

## Matrix example with data from Work et al 2010
data(work)
work_rare <- indiv_rare(as.matrix(work))
### Plot results
plot(x = unlist(lapply(work_rare,length)),unlist(lapply(work_rare,max)),ylim=c(0,max(unlist(work_rare))),xlim=c(0,max(unlist(lapply(work_rare,length)))), ylab = "Species", xlab ="Individuals")
for(i in 1:length(work_rare)){
lines(work_rare[[i]])
}

}
}

