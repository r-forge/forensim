\name{N2error}
\alias{N2error}
\title{
Calculates exact error for maximum allele count for two markers
}
\description{
The maximum allele count principle leads two wrong conclusion for 
two contributors if only a
maximum of one or two alleles is seen. This probability of error is calculated. 
}
\usage{
N2error(dat)
}
\arguments{
  \item{dat}{
 \code{dat} Data frame. First column allele size, remaining columns frequencies
}
}
\value{
The probabilty of error is returned
}
\author{ Thore Egeland \email{Thore.Egeland@medisin.uio.no}}

\examples{
#Example based on 15 markers of Tu data
library(forensim)
data(Tu)
N2error(Tu)
}
