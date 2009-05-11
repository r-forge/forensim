\encoding{UTF-8}
\name{nball}
\alias{nball}
\title{ Number of alleles in a mixture}
\description{
  \code{nball} gives the number of alleles of a  \code{simumix} object. 
}
\usage{
nball(mix, byloc = FALSE)
}

\arguments{
  \item{mix}{ a \code{simumix} object}
  \item{byloc}{a logical indicating whether the number of alleles must be calculated by locus  or for all loci (default)}
}



\author{Hinda Haned  <haned@biomserv.univ-lyon1.fr>}


\seealso{  \code{\link{simumix}} }
\examples{
data(strusa)
#simulating 100 genotypes with allele frequencies from the   African American population
gaa<-simugeno(strusa,n=c(100,0,0))
#simulating a 4-person mixture
maa4<-simumix(gaa,ncontri=c(4,0,0))
nball(maa4,byloc=TRUE)
}

\keyword{ misc }
