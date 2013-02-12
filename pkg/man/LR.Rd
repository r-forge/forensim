\name{LR}
\alias{LR}
\title{Likelihood ratio for DNA evidence interpretation (2): a sophistacted version of function LR()}


\description{LR  Allows the calculation of likelihood ratios for a piece of DNA evidence, for any number of replicates, any number of contributors, and when drop-in and drop-out are possible.
}
\usage{
LR(Repliste, Tp, Td, Vp, Vd, xp, xd, theta, prDHet, prDHom, prC, freq)
}
\arguments{

  \item{Repliste}{vector of alleles present at a given locus for any number of replicates. If there are two replicates, showing alleles 12,13, and 14 respectively, then \code{Repliste} should be given as c(12,13,0,14), where the 0 is used as a separator. An empty replicate is simply 0. For example, replicates (12,13) and and one empty replicate must be given as: c(12,14,0,0). }
  
   
  \item{Tp}{ vector of genotypes for the known contributors under Hp. Genotype 12/17 should be given as a vector c(12,17)
 and genotypes 12/17,14/16, should be given as a unique vector: c(12,17,14,16).}
  \item{Td}{vector of genotypes for the known contributors under  Hd. Should be in the same format as Tp. If there are no known contributors under Hd, then set Td to 0.}
  
  \item{Vp}{vector of genotypes for the known non-contributors (see References section) under Hp. See \code{Tp} for format. 
}
  \item{Vd}{vector of genotypes for the known non-contributors (see References section) under Hd. Should be in the same format than Vp, if empty, set to 0.
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




  \item{prC}{probability of drop-in applied per locus}
  \item{freq}{
  vector of the corresponding allele frequencies of the analysed locus in the target population}
}


\references{

Gill, P.; Kirkham, A. & Curran, J. LoComatioN: A software tool for the analysis of low copy number DNA profiles \emph{Forensic Science International}, 2007, 166(2-3), 128-138

Curran, J. M.; Gill, P. & Bill, M. R. Interpretation of repeat measurement DNA evidence allowing for multiple contributors and population substructure \emph{Forensic Science International}, 2005, 148, 47-53
}


\author{Hinda Haned \email{h.haned@nfi.minvenj.nl}}




\seealso{
 \code{\link{LRmixTK}}

}


\examples{
#load allele frequencies
library(forensim)
data(ngm)
#create vector of allele frequencies
d10<-ngm$tab$D10
# heterozygote dropout probability (resp. homozygote) is set to 0.2 for all
# contributors (0.04 for homozygotes) 
LR(Repliste=c(12,13,14),Tp=c(12,13),Td=0,Vp=0,Vd=0,xd=2,xp=1,theta=0,prDHet=c(0.2,0.2),
prDHom=c(0.04,0.04),prC=0,freq=d10)


}
\keyword{htest}

