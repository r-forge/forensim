\encoding{UTF-8}
\name{Cmn}
\alias{Cmn}
\title{ The number of all possible combinations of m elements among n with repetitions  }
\description{
 The number of all possible combinations of m elements among n with repetitions.
}
\usage{
Cmn(m, n)
}

\arguments{
  \item{m}{ the  \code{m} elements to combine among \code{n} }
  \item{n}{ the  \code{n} elements from which to combine \code{m} elements with repetitions }
}
\details{
 There are (n+m-1)!/(m!(n-1)!) ways to combine m elements among n with repetitions.}
\note{Cmn was implemented as an auxiliary function for the \code{\link{dataL}} function which computes the
likelihood of a mixture alleles conditional on the number of contributors.}
 
\author{Hinda Haned  <haned@biomserv.univ-lyon1.fr>}


\seealso{ \code{\link{comb}} for all possible combinations of m elements among n with repetitions}
\examples{
Cmn(2,3)
comb(2,3)
}

\keyword{ models }

