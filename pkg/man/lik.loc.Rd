\encoding{UTF-8}
\name{lik.loc}
\alias{lik.loc}
\title{Likelihood  per locus of the observed alleles in a DNA mixture conditional on the number of contributors
to the mixture}
\description{
The \code{lik.loc} function computes the likelihood of the observed data in a forensic DNA mixture, for each of the loci involved, conditional on the number of contributors to 
the mixture. 
}
\usage{
lik.loc(x = 1, mix, freq, refpop = NULL, theta = NULL, loc=NULL) 
}

\arguments{
  \item{x}{the number of  contributors to the DNA  mixture }
  \item{mix}{ a \code{simumix} object which contains the mixture to be analyzed }
  \item{freq}{ a \code{tabfreq} object from which to extract the allele frequencies }
  \item{refpop}{ a factor giving the reference population in \code{tabfreq} from which to 
  extract the allele frequencies
  }
  \item{theta}{a float from  [0,1[ giving Wright's  Fst coefficien. \code{theta}
  acounts for population subdivision while computing the likelihood of the data. }
  \item{loc}{the loci for which the likelihood shall be computed. Default (set to NULL) corresponds to all loci. }
}


\details{
\code{lik.loc} computes the likelihood per locus of the observed alleles.
This function implements the general formula for the interpretation
of DNA mixtures in case of subdivided populations (Curran et al, 1999), in the particular case where all contributors 
are unknown and belong to the same subpopulation.\cr
The Fst coefficient given in the \code{theta} argument allows accounting for population subdivision when all
contributors belong to the same subpopulation.\cr
}

\value{
The function \code{lik.loc} returns a vector, of length the number of loci in \code{loc}, 
giving the likelihood of the data for each locus.
}

\references{
Haned H, Pontier D, Lobry JR, Pene L, Dufour AB.
Estimating the number of contributors to forensic DNA mixtures: Does maximum likelihood 
perform better than maximum allele count? Submitted, 2009.\cr

 Curran JM, Triggs CM, Buckleton J, Weir BS.
Interpreting DNA Mixtures in Structured Populations. \emph{ J Forensic Sci} 1999;44(5): 987-995
}

\author{Hinda Haned  <haned@biomserv.univ-lyon1.fr>}




\seealso{  \code{\link{lik}} for the overall loci likelihood, \code{\link{likestim}} and
\code{\link{likestim.loc}} for the estimation of the number of contributors to a DNA mixture through
likelihood maximization}


\examples{
data(strusa)
#simulation of 1000 genotypes from the Caucasian allele frequencies
gen<-simugeno(strusa,n=c(0,100,0))

#4-person mixture
mix4 <- simumix(gen,ncontri=c(0,4,0))
lik.loc(x=2,mix4, strusa, refpop="Cauc")
lik.loc(x=2,mix4, strusa, refpop="Afri")
#You may also want to try:
#likestim(mix4,strusa,refpop="Cauc")

}


\keyword{htest}