\name{rhs}
\alias{rhs}
\title{the right hand side of eqn 4 in Colwell 2012}
\usage{
  rhs(vec, n, m, kclass, freqk)
}
\arguments{
  \item{vec}{the species abundance vector}

  \item{n}{the total number of individuals in the sample}

  \item{m}{the number of individuals to interpolate S for.}

  \item{kclass}{the number of individuals in a given
  frequency class}

  \item{freqk}{a vector of frequency classes}
}
\value{
  the full rhs of eqn 4 as a scalar
}
\description{
  calculate the right hand side of eqn 4, works by only
  calculating values where frequency classes exist and
  summing them
}
\details{
  all these values can be calculated from the vector of
  species in the initial function call but I pass them for
  transparency's sake
}

