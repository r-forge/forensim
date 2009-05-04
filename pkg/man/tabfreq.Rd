\encoding{UTF-8}
\name{tabfreq}
\docType{class}
\alias{tabfreq-class}
\alias{names,tabfreq-method}
\alias{print,tabfreq-method}
\alias{show,tabfreq-method}

\title{forensim class for population allele frequencies}
\description{	The S4 \code{tabfreq} class is used to store allele frequencies, from either one or several populations.
}

\section{Slots}{
	 \describe{
    \item{\code{tab}:}{   a list giving allele frequencies for each locus. If there are several populations,
	\code{tab}	gives allele frequencies  in each population  }
     \item{\code{which.loc}:}{ character vector giving the names of the loci}
     \item{\code{pop.names}:}{  factor of populations names (optional)}
	
  }}




\section{Methods}{
  \describe{
    \item{names}{\code{signature(x = "tabfreq")}}: gives the names of the attributes of a tabfreq object
    \item{show}{\code{signature(object = "tabfreq")}}:  shows a tabfreq object
	\item{print}{\code{signature(object="tabfreq")}}: prints a tabfreq object
	 
	 }
}




\seealso{
	 \code{\link{as.tabfreq}}, \code{\link{is.tabfreq}} and \code{\link{simugeno}} for genotypes simulation from allele frequencies stored in a \code{tabfreq} object} 



\author{ Hinda Haned  \email{haned@biomserv.univ-lyon1.fr} }


\examples{
showClass("tabfreq")
}
\keyword{classes}
\keyword{manip}
\keyword{datagen}
