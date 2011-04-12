\encoding{UTF-8}
\name{simufreqD}
\alias{simufreqD}

\title{ Function to simulate allele frequencies for independent loci from  a Dirichlet model}
\description{
 The \code{simufreqD} function simulate single population allele frequencies for independent loci.
 Allele frequencies are generated as random deviates from a Dirichlet distribution, whose parameters control
 the mean and the  variance of the  simulated allele frequencies. 
}
\usage{
simufreqD(nloc = 1, nal = 2, alpha = 1)
}

\arguments{
  \item{nloc}{ the number of loci to simulate}
  \item{nal}{ the numbers of alleles per locus. Either an integer, if the loci have the same number of 
  alleles, or an integer  vector, if the number of alleles differ between loci}
  \item{alpha}{ the parameter used to simulate allele frequencies from the Dirichlet distribution. If the 
  \code{nloc} loci have the same allele number, 
  \code{alpha} can either be the same for all alleles (default is one: uniform distribution), in this case \code{alpha} is an integer, or 
  \code{alpha} can  be  different  between alleles at a given locus, in this case, \code{alpha} is a matrix of dimension \code{nal} x \code{nloc}. 
  
  When  the number of alleles differ between loci, \code{alpha} can either be the same or differ between alleles  at a given locus.
  In the first  case \code{alpha} is a vector of length \code{nloc}, in the second case,  \code{alpha} is a matrix of dimensions
 \code{nal} x \code{nloc} where NAs are introduced for alleles not seen at a given locus.
 }
}


\details{
Allele frequencies for independent loci  are simulated using a Dirichlet distribution with parameter 
\code{alpha}. At  a given locus  L with n alleles,  the  allele frequencies are modeled as a  vector of random 
variables
p=(p1, ..., pn), following a Dirichlet distribution with parameters:\cr
alpha = (alpha1, ..., alphan) where p1+...+pn=1  and  alpha1,..., alphan  > 0. 

}

\value{
  
  A matrix containing the simulated allele frequencies. The data is presented in the format of the Journal of
  Forensic Sciences for genetic data: allele names are given in the first column, and frequencies for a given allele are read in rows for the
  different markers in columns. When an allele is not observed for a given locus, the value is coded  NA (instead of  "-" in the 
  original format).
}
\references{ 
Johnson NL, Kotz S, Balakrishnan N. Continuous Univariate Distributions, vol 2. John Wiley & Sons, 1995.\cr

Wright S. The genetical structure of populations. Ann Eugen 1951;15:323-354. 
}


\author{Hinda Haned \email{h.haned@nfi.minjus.nl}}
\note{ 
 The code used here  for the generation of random Dirichlet deviates was previously implemented in the 
 gtools library.}
\seealso{ \code{\link{simupopD}}}

\examples{
#simulate alleles frequencies for 5 markers with respectively 2, 3, 4, 5, and 6 alleles

simufreqD(nloc=5,na=c(2,3,4,5,6) , alpha=1)
}


\keyword{datagen }

