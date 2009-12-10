\name{N2Exact}
\alias{N2Exact}

\title{
Calculates exact allele distribution for 2 contributors
}
\description{
The distribution of N, the number of alleles
showing is calculated exactly assuming 2 contributors.
Theta-correction is not implemented.
The function may be used to check accuracy of simulations and
indicate required number of simulations for one example.
}
\usage{
N2Exact(p)
}
\arguments{
  \item{p}{
vector of allele frequencies. Must sum to 1. Default: for uniformly distrubted alleles.
}
}

\value{
Returns(P(N=i) for i=1,2,3,4
}

\author{ Thore Egeland \email{Thore.Egeland@medisin.uio.no}}


%\note{
%The error involved for maximum allele count roule can be calculated
%using this function since the probabilyt that the maximum is <=2 equals a product 
%}
\examples{
#Distribution for a marker with 20 alles of equal frequency
N2Exact(p=rep(0.05,20))
}
