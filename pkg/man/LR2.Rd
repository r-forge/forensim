\name{LR2}
\alias{LR2}
\title{ Likelihood ratio for DNA evidence interpretation (2): a sophistacted version of function LR()}


\description{ LR2  Allows the calculation of likelihood ratios for a piece of DNA evidence, for any number of replicates, any number of contributors, and when drop-in and drop-out are possible.
}
\usage{
LR2(Repliste, Tp, Td, Vp, Vd, xp, xd, theta, prDHet, prDHom, prC, freq)
}
\arguments{
  \item{Repliste}{list of the alleles present at a given locus for any number of replicates. If there are two replicates, showing alleles 12,13, and 14 respectively, then \code{Repliste} should be given as  list(c(12,13),14).}
  
  
  \item{Tp}{ vector of genotypes for the known contributors under a Hp. Should be given in the genetics format, i.e., genotype 1217 should be given as '12/17' (character), and genotypes 1217,1416, should be given as a vector: c('12/17','14/16').}
  \item{Td}{vector of genotypes for the known contributors under a Hd. Should be in the same format as Tp. If there are no known contributors under Hd, then set Td to NULL.}
  
  \item{Vp}{vector of genotypes for the known non-contributors (see References section) under a Hp. Should be given in the genetics format, i.e., genotype 1217 should be given as '12/17' (character), and genotypes 1217,1416, should be given as a vector: c('12/17','14/16') 
}
  \item{Vd}{vector of genotypes for the known non-contributors (see References section) under a Hd. Should be in the same format than Vp, if empty, set to NULL.
}
  \item{xp}{Number of unknown individuals under Hd. Set to 0 if there are no unknown contributors.}

  \item{xd}{Number of unknown individuals under Hd. Set to 0 if there are no unknown contributors.}
  \item{theta}{thete correction, value must be taken in [0,1)}


  \item{prDHet}{
probability of dropout for heterozygotes. It is possible to assign different values per contributor. In this case, \code{prDHet} must be a vector, of length the number of contributors in T + x, and the probabilities must be given in this order. if the probability of dropout for T is d1, and for the unknown is d2, then \code{prDHet}=(d1,d2). In case T is not a heterozygote, the given vector must still be of length length(T) +x, but the given value for T does not matter, because it won't be used,
the value in prDHom is used instead. This is a bit ad hoc and an improvement is currently under development. }
  \item{prDHom}{
probability of dropout for homozygotes. See description ofr argument \code{PrDHom}.
}
  \item{prC}{
probability of drop-in}
  \item{freq}{
  vector of the corresponding allele frequencies of the analysed locus in the target population}
}


\references{

Gill, P.; Kirkham, A. & Curran, J. LoComatioN: A software tool for the analysis of low copy number DNA profiles \emph{Forensic Science International}, 2007, 166(2-3), 128-138

Curran, J. M.; Gill, P. & Bill, M. R. Interpretation of repeat measurement DNA evidence allowing for multiple contributors and population substructure \emph{Forensic Science International}, 2005, 148, 47-53
}


\author{Hinda Haned \email{h.haned@nfi.minjus.nl}}




\seealso{
 \code{\link{LRmixTK}}

}


\examples{
#load allele frequencies
library(forensim)
data(ngm)
#create vector of allele frequencies
d10<-ngm$tab$D10
#single  dropout parameter for all contributors
LR2(list(12,13,14),Tp='12/13',Td=NULL,Vp=NULL,Vd=NULL,xd=2,xp=1,theta=0,prDHet=0.2,prDHom=0.1,prC=0,freq=d10)

#different dropout probabilities for contributor 1 in T and the unknown contributor: 0.2 for the first one (heterozygote, corresponding value is 0 in the homozygote vector), 0.7 for the second one if it is a heterozygote and 0.3 if its is a homozygote.

LR2(list(12,13,14),Tp='12/13',Td=NULL,Vp=NULL,Vd=NULL,xd=2,xp=1,theta=0,prDHet=c(0.2,0.70),prDHom=c(0,0.3),prC=0,freq=d10)
}
\keyword{htest}

