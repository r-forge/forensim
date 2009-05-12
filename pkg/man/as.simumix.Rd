\encoding{UTF-8}
\name{simumix constructor}
\alias{simumix-methods}
\alias{simumix}
\alias{as.simumix}
\alias{is.simumix}
\title{simumix constructor}
\description{Constructor for \linkS4class{simumix} objects.\cr
  The function \code{simumix} creates a \linkS4class{simumix} object from
  a \linkS4class{tabfreq} object.\cr

  The function \code{as.simumix} is an alias for \code{simumix} function.\cr
 
  \code{is.simumix} tests if an object is a valid simumix object.\cr

  Note: to get the manpage about \linkS4class{simumix}, please type
  'class ? simumix'.
}

\usage{
simumix(tab,which.loc=NULL,ncontri=1)
as.simumix(tab,which.loc=NULL,ncontri=1)
is.simumix(x)
}

\arguments{
  \item{tab}{a simugeno object created with constructor simugeno}
  \item{which.loc}{a character vector giving the chosen loci for the genotypes simulation. The default is set to NULL, 
  which corresponds to all the loci of the \code{simugeno} object given in argument }
  \item{ncontri}{integer vector giving the number of individuals. If there are several populations, 
  \code{ncontri} gives the numbers of individuals to simulate per population. Default is one.}
  \item{x}{an object}}


\value{
  For \code{simumix} and \code{as.simumix}, a simumix object. For \code{is.simumix}, a logical.
}

\seealso{
  \code{"\linkS4class{simumix}"}, \code{\link{simugeno}} for creating a simugeno object.
}


\author{Hinda Haned \email{haned@biomserv.univ-lyon1.fr}}
\examples{
data(Tu)
tab<-simugeno(tabfreq(Tu),n=1200)
#simulation of a 3-person mixture characterized with markers FGA, TH01 and TPOX 
simumix(tab,which.loc=c('FGA','TH01', 'TPOX') , n =3)

}
\keyword{manip}
\keyword{datagen}