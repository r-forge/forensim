\encoding{UTF-8}
\name{findfreq}
\alias{findfreq}

\title{Finds the frequencies of the alleles of a mixture from a tabfreq object }
\description{
 The \code{findfreq} function finds the frequencies of the alleles of a mixture stored in a simumix object, form a given tabfreq object. 
 If the tabfreq object contains multiple populations, a reference population from which to extract the  
 frequencies must be specified.
}
\usage{
findfreq(mix, freq, refpop = NULL)
}
\arguments{
  \item{mix}{ a \code{simumix} object }
  \item{freq}{ a \code{tabfreq} object from which to extract the alleles frequencies of the mixture }
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

