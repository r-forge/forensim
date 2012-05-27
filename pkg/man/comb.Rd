\encoding{UTF-8}
\name{comb}
\alias{comb}
\title{ Generate all possible combinations of  m elements among n with repetitions  }
\description{
  Generate all possible combinations of m elements among n with repetitions. 
}
\usage{
comb(m, n)
}

\arguments{
  \item{m}{ the number of elements to combine}
  \item{n}{ the number of elements from which to combine the \code{m} elements}
}
\details{
 There are (n+m-1)!/(m!(n-1)!) ways to combine m elements among n with repetitions, \code{combn} generates 
 all these possible combinations.
 
}
\value{
A matrix of  (n+m-1)!/(m!(n-1)!) rows, and n columns, each row is  a possible combination of  m elements among n .
}
  

\author{Hinda Haned \email{h.haned@nfi.minvenj.nl}}

\seealso{  \code{\link{Cmn}} for the calculation of the number of all possible combinations of m elements among n with repetitions  }
\examples{
#combine 2 objects among 3 with repetitions
Cmn(2,3)
comb(2,3)
}
\keyword{ models }
