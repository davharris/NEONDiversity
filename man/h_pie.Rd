\name{h_pie}
\alias{h_pie}
\title{Calculate Hurlbert's PIE}
\usage{
  h_pie(sp_list)
}
\arguments{
  \item{sp_list}{a vector of species abundances, or a site
  (rows) x species (columns) matrix}
}
\value{
  a scalar (for a vector of input) or vector (for matrix
  input) of PIE values
}
\description{
  Calculate Hurlbert's PIE (Probability of an Interspecific
  Encounter)
}
\details{
  This will calculate Hurlbert's PIE after Hurlbert 1971.
  The equation used is \eqn{PIE = \begin{bmatrix}
  \frac{N}{N-1} \end{bmatrix} \begin{bmatrix} 1 -
  \sum_{i}^{N}\begin{pmatrix} \frac{N_i}{N} \end{pmatrix}^2
  \end{bmatrix}}{ PIE = \begin{bmatrix} \frac{N}{N-1} }
}
\examples{
\dontrun{
## Generate some data
fake_dat <- matrix(rpois(100,100),ncol=10,nrow=10)

### Return 10 values for PIE

h_pie(fake_dat)

}
}

