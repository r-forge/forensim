\name{PE}
\Rdversion{1.1}
\alias{PE}

\title{
Exclusion probability}
\description{
Computes the exclusion probability of a mixture stored in a \code{\linkS4class{simumix}} object
}
\usage{
PE(mix, freq, refpop = NULL, theta = 0, byloc = FALSE)
}

\arguments{
\item{mix}{ a \code{simumix} object }
\item{freq}{ a \code{tabfreq} object giving the allele frequencies from which to compute the
exclusion probability}
\item{refpop}{character giving the reference population, used only if \code{freq} contains 
allele frequencies for multiple populations
}
\item{theta}{a float from [0,1[ giving Wright's  Fst coefficient. \code{theta}
 accounts for population subdivision while computing the likelihood of the data. }
\item{byloc}{logical, if TRUE, than the exclusion probability is computed per locus, if FALSE (default),
the calculations are done for all loci simultaneously
}

}
\details{
\code{PE} gives the exclusion probability  at  a locus, or at several loci when conditions for Hardy Weinberg are
met. If this condition is not met in the population,
than a value for \code{theta} must be supplied to take into account dependencies
between alleles. The formula of the exclusion probability that allows taking into account departure
 from Hardy  Weinberg proportions due to population subdivision  was provided by Bruce Weir, please see the
 references section.
}

\references{
Clayton T, Buckleton JS. Mixtures. In: Buckleton JS, Triggs CM, Walsh SJ, editors.
Forensic DNA Interpretation. CRC Press 2005;217-74

}

\author{Hinda Haned  <haned@biomserv.univ-lyon1.fr>}




\examples{
data(strusa)
geno1<-simugeno(strusa,n=c(0,0,100))
mix2 <-simumix(geno1,ncontri=c(0,0,2))
PE(mix2,strusa,"Hisp",byloc=TRUE)
}
\keyword{htest}
