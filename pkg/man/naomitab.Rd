\encoding{UTF-8}
\name{naomitab}
\alias{naomitab}

\title{ Handling of missing values in a data frame }
\description{
  \code{naomitab} handles missing values (NA) in a data frame: it returns  a list  of the columns where NA's have
  been  removed.
}
\usage{
naomitab(tab)
}

\arguments{
  \item{tab}{ a data frame}
}

\value{
Returns a list of  length the number of columns in \code{tab}
where each component is a column of \code{tab}, and the values are the corresponding rows where NA's have
been removed.
  
}
\note{
 This function was designed to handle missing values in data frames in the format of the Journal of Forensic
 Sciences for population genetic data:  allele names are given in the first column, and frequencies for a 
 given allele   are read in rows for  different loci. 
  When a given allele is not observed, the value is  encoded NA (originally encoded  "-" in the journal).
}
\author{Hinda Haned  <haned@biomserv.univ-lyon1.fr>}

\seealso{ \code{\link{tabfreq}} }

\examples{
data(Tu)
naomitab(Tu)
}

\keyword{ manip }
