\encoding{UTF-8}
\name{simugeno}
\alias{simugeno-class}
\alias{names,simugeno-method}
\alias{print,simugeno-method}
\alias{show,simugeno-method}

\title{forensim class for simluated genotypes}
\description{The S4 \code{simugeno} class is used to store existing or simulated genotypes.}


\section{Slots}{
	 \describe{
    \item{\code{tab.freq}:}{  a list giving allele frequencies for each locus.
	If there are several populations,
	\code{tab.freq}	gives allele frequencies  in each population }
    \item{\code{nind}:}{ integer vector giving the number of individuals. If there are several
	populations, \code{nind} gives the numbers of  individuals per population}
    \item{\code{pop.names}:}{ factor of populations names}
    \item{\code{popind}:}{ factor giving the population of each individual   }
    \item{\code{which.loc}:}{character vector giving the locus names}
    \item{\code{tab.geno}:}{ matrix giving the genotypes (in rows) for each locus (in columns). 
	The genotype of a homozygous individual carrying the allele "12" is coded "12/12". A heterozygous 
			individual carrying alleles "12" and "13" is coded "12/13" or "13/12".}
    \item{\code{indID}:}{ character vector giving the individuals ID}
  }
}

\section{Methods}{
  \describe{
    \item{names}{\code{signature(x = "simugeno")}: gives the names of the attributes of 
	a simugeno object }
    \item{show}{\code{signature(object = "simugeno")}: shows a simugeno object }
	\item{print}{\code{signature(object = "simugeno")}: prints a simugeno object }
	 }
}



\seealso{\code{\link{as.simugeno}} for the simugeno class constructor,
  \code{\link{is.simugeno}}, \code{\link{simumix}} and 
\code{\link{tabfreq}}}  


\author{ Hinda Haned  \email{contact@hindahaned.info} }



\examples{
showClass("simugeno")

}
\keyword{classes}
\keyword{manip}
\keyword{datagen}
