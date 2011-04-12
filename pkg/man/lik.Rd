\encoding{UTF-8}
\name{lik}
\alias{lik}
\title{ Likelihood of the observed alleles at different loci in a DNA mixture conditional on the number of 
contributors to the mixture}
\description{
The \code{lik} function computes the likelihood of the observed alleles in a forensic DNA mixture, for a set 
of loci, conditional on the number of contributors to the mixture. The overall likelihood is computed as the 
product of loci likelihoods.
}

\usage{
lik(x = 1, mix, freq, refpop = NULL, theta = NULL, loc=NULL)
}

\arguments{
  \item{x}{the number of contributors to the DNA mixture, default is 1 }
  \item{mix}{ a \code{simumix} object which contains the mixture to be analyzed  }
  \item{freq}{ a \code{tabfreq} object from which to extract the allele frequencies }
  \item{refpop}{ a factor giving the reference population in \code{tabfreq} from which to 
  extract the allele frequencies. This argument is used only if \code{freq} contains allele frequencies for multiple
  populations, otherwise it is by default set to NULL
  }
  \item{theta}{a float from  [0,1[ giving Wright's  Fst coefficient.
  \code{theta} accounts for population subdivision while computing the likelihood of the data }
  \item{loc}{ loci for which the overall likelihood shall be computed.  Default (NULL) corresponds to all loci}
}


\details{
\code{lik} computes the likelihood of the alleles observed at all loci conditional on the number of contributors. 
  This function implements the general formula for the interpretation of DNA mixtures
 in case of population subdivision (Curran et al, 1999), in the particular case where all contributors are unknown
 and belong to the same subpopulation.\cr
  The likelihood for multiple loci is computed as the product of loci likelihoods.
}




\references{ 
Haned H, Pene L, Lobry JR, Dufour AB, Pontier D.
Estimating the number of contributors to forensic DNA mixtures: Does maximum likelihood 
perform better than maximum allele count? \emph{ J Forensic Sci}, accepted 2010.\cr
Curran JM, Triggs CM, Buckleton J, Weir BS. Interpreting DNA Mixtures in Structured Populations. 
\emph{ J Forensic Sci} 1999;44(5): 987-995


}
\author{Hinda Haned  <h.haned@nfi.minjus.nl>}



\seealso{  \code{\link{lik.loc}} for the likelihood per locus, \code{\link{likestim}} and
\code{\link{likestim.loc}} for the estimation of the number of contributors to a DNA mixture through
likelihood maximization}
\examples{
data(strusa)
#simulation of 1000 genotypes from the African American allele frequencies
gen<-simugeno(strusa,n=c(1000,0,0))
#3-person mixture
mix3<-simumix(gen,ncontri=c(3,0,0))
sapply(1:3, function(i) lik(x=i,mix3, strusa, refpop="Afri"))

}


\keyword{htest}