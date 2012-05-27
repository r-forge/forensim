\encoding{UTF-8}
\name{dataL}
\alias{dataL}
\title{Generic formula of the likelihood of the observed alleles in a mixture conditional on the number 
of contributors for a specific locus}
\description{
 The function \code{dataL} gives the likelihood of a set of alleles observed at a specific locus conditional 
 on the number of contributors that gave  these alleles.  Calculation is based upon the frequencies
of the observed alleles.


}
\usage{
dataL(x = 1, p, theta = 0)
}

\arguments{
  \item{x}{  an integer giving the number of contributors }
  \item{p}{ a  numeric vector giving the frequencies of the observed alleles in the mixture}
  \item{theta}{ a float  in [0,1[. \code{theta} is equivalent to Wright's Fst. In case of population
  subdivision, it allows a correction of the allele frequencies in the subpopulation of interest}}


\references{ 

Haned H, Pene L, Lobry JR, Dufour AB, Pontier D.
Estimating the number of contributors to forensic DNA mixtures: Does maximum likelihood 
perform better than maximum allele count? \emph{ J Forensic Sci}, accepted 2010.\cr

Curran JM, Triggs CM, Buckleton J,  Weir BS.
Interpreting DNA Mixtures in Structured Populations. \emph{ J Forensic Sci} 1999;44(5): 987-995}

\note{\code{dataL} function has several similarities with the \code{Pevid.gen} function 
 of the \emph{forensic} package which computes the probability of the DNA evidence, \code{dataL}
implements a particular case of this probability. Please see \url{http://cran.r-project.org/web/packages/forensic/}}
\author{Hinda Haned \email{h.haned@nfi.minvenj.nl}}

\seealso{\code{\link{lik.loc}} and  \code{\link{lik}} for calculating the likelihood of a given simumix object }
\examples{
#likelihood of observing two alleles at frequencies 0.1 and 0.01 when the number of 
#contributors is 2, in two cases:  theta=0 and theta=0.03
dataL(x=2,p=c(0.1,0.01), theta=0)
dataL(x=2,p=c(0.1,0.01), theta=0.03)
}

\keyword{htest}

