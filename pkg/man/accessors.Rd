\encoding{UTF-8}
\docType{methods}
\name{Accessors}
\alias{$,simugeno-method}
\alias{$,simumix-method}
\alias{$,tabfreq-method}
\alias{$<-,simugeno-method}
\alias{$<-,simumix-method}
\alias{$<-,tabfreq-method}

\title{ Accessors for forensim objects}
\description{
Accessors for forensim objects:  \linkS4class{simugeno}, \linkS4class{simumix} and \linkS4class{tabfreq}.
"\$" and "\$<-" are used to access the slots of an object, they are equivalent 
to "@" and "@<-".\cr
}

\value{
  A \linkS4class{simugeno}, a \linkS4class{simumix} or a \linkS4class{tabfreq} object.
}

\author{ Hinda Haned \email{haned@biomserv.univ-lyon1.fr} }



\examples{
data(strusa)
class(strusa)

strusa@pop.names
#equivalent 
strusa$pop.names
}
\keyword{manip}