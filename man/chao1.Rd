\name{chao1}
\alias{chao1}
\title{Chao1 estimator}
\usage{
  chao1(spp = NULL, f1 = NULL, f2 = NULL)
}
\arguments{
  \item{spp}{a vector or a site x species matrix}

  \item{f1}{a scalar or vector of the frequency of
  singletons see details}

  \item{f2}{a scalar or vector of the frequency of
  doubletons see details}
}
\value{
  a vector of estimates unobserved species which can be
  added to S for an estimate of the total number of species
}
\description{
  Calculates the Chao1 estimator based on Colwell 2012 and
  Chao 1984
}
\details{
  The function can take either a site by species matrix or
  alternatively the parameters f1, f2 and S.  Those are
  calculated from the matrix / vector, but don't need to be
  if that is all the information you have
}
\examples{
\dontrun{
data(costaRicaTrees)
chao1(new_growth)

}
}

