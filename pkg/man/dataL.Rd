\encoding{UTF-8}
\name{dataL}
\alias{dataL}
\title{Generic formula of the likelihood of the observed alleles in a mixture conditional on the number of contributors for a specific locus}
\description{
 The function \code{dataL} gives the likelihood of a set of alleles observed at a specific locus, given as allele
 frequencies, conditionnal on the number of contributors that gave  the observed alleles.
}
\usage{
dataL(x = 1, p, theta = 0)
}

\arguments{
  \item{x}{  the uknown number of contributors }
  \item{p}{ a vector giving the allele frequencies of the observed alleles in the mixture}
  \item{theta}{ an integer  in [0,1[. \code{theta} is the correction factor for   population subdivision
}}


\references{ Curran JM, Triggs CM, Buckleton J,  Weir BS.
Interpreting DNA Mixtures in Structured Populations. \emph{ J Forensic Sci} 1999;44(5): 987-995}
\note{\code{dataL} function has several similarities with the \code{Pevid.gen} function 
 of the \emph{forensic} package which gives the match probability of the evidentiary DNA profile with a
 suspect's profile, \code{dataL}  implements a particular case of this probability. Please see \url{http://cran.r-project.org/web/packages/forensic/}}
\author{Hinda Haned \email{haned@biomserv.univ-lyon1.fr}}

\seealso{\code{\link{lik.loc}} and  \code{\link{lik}} for calculating the likelihood of a given simumix object }
\examples{
#likelihood of observing two at frequencies 0.1 and 0.01 when the number of contributors is 2
#in two cases: Wright's Fst=0 and Fst=0.03
dataL(x=2,p=c(0.1,0.01), theta=0)
dataL(x=2,p=c(0.1,0.01), theta=0.01)
}

\keyword{htest}

