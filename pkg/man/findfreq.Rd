\encoding{UTF-8}
\name{findfreq}
\alias{findfreq}

\title{Find the frequencies of a mixture's alleles from a tabfreq object }
\description{
 Find the frequencies of a  mixture' alleles  form a given tabfreq object. If multiple populations are stored
 in the tabfreq object, a reference population, from which to extract the a frequencies must be specified.
}
\usage{
findfreq(mix, freq, refpop = NULL)
}
\arguments{
  \item{mix}{ a \code{simumix} object }
  \item{freq}{ a \code{tabfreq} from which to extract the mixture allele frequencies }
  \item{refpop}{a factor giving the reference population in \code{tabfreq} from which to extract the 
  allele frequencies }
}





\value{
A list giving the allele frequencies for each locus.
}
\author{Hinda Haned  <haned@biomserv.univ-lyon1.fr>}

\seealso{\code{\link{simumix}}}
\examples{
data(strusa)
s2<-simumix(simugeno(strusa,n=c(0,2000,0)),ncontri=c(0,2,0))
findfreq(s2,strusa,refpop="Cauc")
}

\keyword{ misc }

