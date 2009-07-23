\encoding{UTF-8}
\name{findmax}
\alias{findmax}

\title{ Function to find the maximum of a vector and its position}
\description{
  The \code{findmax} function finds the maximum of a vector and its position.
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
A matrix of two columns:\cr
\code{max} the position of the maximum in \code{vec} \cr
\code{maxval} the  maximum 
}

\note{\code{findmax} is an auxiliary function for the \code{\link{dataL}} function, 
used to compute the
likelihood of the observed alleles in a mixed DNA stain given the number of contributors.  }

\author{Hinda Haned  <haned@biomserv.univ-lyon1.fr>}



\examples{

findmax(1:10)
}

\keyword{misc}
