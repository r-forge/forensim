\encoding{UTF-8}
\name{simugeno constructor}
\alias{simugeno-methods}
\alias{simugeno}
\alias{as.simugeno}
\alias{is.simugeno}
\title{simugeno constructor}
\description{Constructor for \linkS4class{simugeno} objects.\cr
  The function \code{simugeno} creates a \linkS4class{simugeno} object from
  a \linkS4class{tabfreq} object.\cr

  The function \code{as.simugeno} is an alias for \code{simugeno} function.\cr
 
  \code{is.simugeno} tests if an object is a valid simugeno object.\cr

  Note: to get the manpage about \linkS4class{simugeno}, please type
  'class ? simugeno'.
}

\usage{
simugeno(tab,which.loc=NULL,n=1)
as.simugeno(tab,which.loc=NULL,n=1)
is.simugeno(x)
}

\arguments{
  \item{tab}{a tabfreq object created with constructor \code{tabfreq}}
  
  \item{which.loc}{a character vector giving the chosen loci for the genotypes simulation. The default is set to NULL, 
  which corresponds to all the loci of the \code{tabfreq} object given in argument }
	
  \item{n}{ integer vector giving the number of individuals. If there are several
	populations, \code{n} gives the numbers of  individuals to simulate per population. For a single population, default is 1.}
  \item{x}{an object}
}

\details{
At a given locus, an individual's genotype is simulated by randomly drawing two 
alleles (with replacement) 
at their respective allele frequencies in the target population.  
}
\value{
  For \code{simugeno} and \code{as.simugeno}, a simugeno object. For \code{is.simugeno}, a logical.
}

\seealso{
  \code{"\linkS4class{simugeno}"}, and \code{\link{tabfreq}} for creating a tabfreq object from a data file.
}


\author{Hinda Haned \email{haned@biomserv.univ-lyon1.fr}}
\examples{
data(Tu)
tab<-tabfreq(Tu)
#simulation of 3 individual genotypes for the STR  marker FGA
geno1 <- simugeno(tab,which.loc='FGA', n =1000)
geno1@tab.geno
}
\keyword{manip}
\keyword{datagen}