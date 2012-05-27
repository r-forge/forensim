\name{changepop}
\alias{changepop}

\title{
Function to change population-related information in forensim objects
}
\description{
The \code{changepop} function changes population-related information in  \linkS4class{tabfreq},  \linkS4class{simugeno} and
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
a character vector giving the population names to be changed
}
  \item{newpop}{
a character vector  giving the new population names
}
}

\value{
a \code{\link{forensim}} object where the slots containing population-related information have been modified
}

\author{Hinda Haned \email{h.haned@nfi.minvenj.nl}}



\examples{
data(strveneto)
tab1 <- simugeno(strveneto,n=100)
tab2 <- changepop(tab1,"Veneto","VENE")
tab1$pop.names
tab2$pop.names
}
\keyword{ manip }