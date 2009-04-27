\name{changepop}
\Rdversion{1.1}
\alias{changepop}

\title{
Function to change population related informations in forensim objects.
}
\description{
\code{changepop} changes population related informations in  \linkS4class{tabfreq},  \linkS4class{simugeno} and
\linkS4class{simumix} objects
}
\usage{
changepop(obj, oldpop, newpop)
}

\arguments{
  \item{obj}{
a forensim object, either a \linkS4class{tabfreq},  a \linkS4class{simugeno} or  a \linkS4class{simumix}
object
}
  \item{oldpop}{
character giving the population names to be changed
}
  \item{newpop}{
character giving the new population names
}
}

\value{
a \link{forensim} object where the slots containing population related informations have been modified
}

\author{Hinda Haned \email{haned@biomserv.univ-lyon1.fr}}



\examples{
data(strveneto)
tab1 <- simugeno(strveneto,n=100)
tab2 <- changepop(tab1,"Veneto","VENE")
tab1$pop.names
tab2$pop.names
}
\keyword{ manip }