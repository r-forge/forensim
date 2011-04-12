\encoding{UTF-8}
\name{tabfreq constructor}
\alias{tabfreq-methods}
\alias{tabfreq}
\alias{as.tabfreq}
\alias{is.tabfreq}
\title{tabfreq constructor}

\description{Constructor for \linkS4class{tabfreq} objects.\cr
  The function \code{tabfreq} creates a \linkS4class{tabfreq} object from
   a data frame or a matrix giving allele frequencies for a single population in the Journal of Forensic Sciences (JFS) format for population genetic data. 
   Whene multiple populations are considered, data shall  be given as a list, where each element is either a matrix or a data frame in the JFS format, and the 
   populations names must be specified.\cr

  The function \code{as.tabfreq} is an alias for the \code{tabfreq} function.\cr
 
  \code{is.tabfreq} tests if an object is a valid tabfreq object.\cr

  Note: to get the manpage about \linkS4class{tabfreq}, please type
  'class ? tabfreq'.
}

\usage{
tabfreq(tab,pop.names=NULL)
as.tabfreq(tab,pop.names=NULL)
is.tabfreq(x)
}

\arguments{
  \item{tab}{either a matrix or a data.frame of markers allele frequencies given in the Journal of Forensic Sciences format for population genetic data}
  \item{pop.names}{(optional) a factor giving the populations names. For a single population in \code{tab}, default is set to NULL.}
  \item{x}{an object}
 
}

\value{
  For \code{tabfreq} and \code{as.tabfreq}, a tabfreq object. For \code{is.tabfreq}, a logical.
}

\seealso{
  \code{"\linkS4class{tabfreq}"}, \code{\link{simugeno}} for creating a simugeno object from a tabfreq object.
}


\author{Hinda Haned \email{h.haned@nfi.minjus.nl}}
\examples{
data(Tu)
tabfreq(Tu,pop.names=factor("Tu"))


}
\keyword{manip}
\keyword{datagen}