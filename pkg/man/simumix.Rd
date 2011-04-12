\encoding{UTF-8}
\name{simumix}
\alias{simumix-class}
\alias{names,simumix-method}
\alias{print,simumix-method}
\alias{show,simumix-method}

\title{forensim class for  DNA mixtures}
\description{The S4 \code{simumix} class is used to store DNA mixtures of individual genotypes
along with informations about the individuals poulations and the loci used to simulate the genotypes.}

\section{Slots}{
	 \describe{
    \item{\code{ncontri}:}{ integer vector giving the number of contributors to the DNA mixture. If there are 
	several populations, \code{ncontri} gives the number of contributors per population}
    \item{\code{mix.prof}:}{ matrix giving the contributors genotypes (in rows) for each locus (in columns). 
	The genotype of a homozygous individual carrying the allele "12" is coded "12/12". A heterozygous 
			individual carrying alleles "12" and "13" is coded "12/13" or "13/12".}
    \item{\code{mix.all}:}{list giving the alleles present in the mixture for each locus}
    \item{\code{which.loc}:}{ character vector giving the locus  names}
    \item{\code{popinfo}:}{ factor giving the population of each contributor }
  }
}

\section{Methods}{
  \describe{
    \item{names}{\code{signature(x = "simumix")}: gives the names of the attributes of a simumix object }
    \item{show}{\code{signature(object = "simumix")}: shows a simumix object}
	\item{print}{\code{signature(object = "simumix")}: prints a simumix object }
	 }
}


\seealso{ \code{\linkS4class{simugeno}}, \code{\link{as.simumix}}, \code{\link{is.simumix}}, \code{\link{simugeno}} and  \code{\link{tabfreq}}}  

\author{ Hinda Haned  \email{h.haned@nfi.minjus.nl} }

\examples{

showClass("simumix")
data(strusa)

}
\keyword{classes}
\keyword{manip}
\keyword{datagen}
