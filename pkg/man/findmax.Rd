\encoding{UTF-8}
\name{findmax}
\alias{findmax}

\title{ Find the maximum of a vector and its position}
\description{
  Gives the maximum of a vector and its position.
}
\usage{
findmax(vec)
}

\arguments{
  \item{vec}{ a numeric vector}
}

\details{
\code{findmax} finds the maximum value of a vector and its position. 
}
\value{
  A matrix of dimension 2 x 1, the first row  \code{max} gives  the position of the maximum value, and the second row 
\code{maxval} returns the maximum value appearing in \code{vec}.  
}
\note{findmax is an auxiliary function for the \link{dataL} function, used to compute the
likelihood of a mixture alleles given the number of contributors.  }

\author{Hinda Haned  <haned@biomserv.univ-lyon1.fr>}



\examples{

findmax(1:10)
}

\keyword{misc}
