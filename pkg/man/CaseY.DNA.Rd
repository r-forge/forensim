\name{CaseY.DNA}
\alias{CaseY.DNA}
\title{
Alleles of mixture, known contributor, defendant as well as missing (drop-out) allelels
for US case (CaseY, see http://www.zebfontaine.eclipse.co.uk/djb.htm)
}
\description{
There are four lines for each marker correspinding to alleles
for the (i) the mixture (ii) known contributor (if ther is one) (iii) the defendant and
(iv) drop-out alleles. The names and order of the markers should be as for dataBase.
}
\usage{data(Bates.DNA)}
\format{
  A data frame with 40 observations on the following 6 variables.
  \describe{
    \item{\code{Marker}}{
arker names}
    \item{\code{Source}}{Alleles of \code{Mixture}, \code{Known} contributor and \code{defendant} and missing (drop-out) alleles}
    \item{\code{A1}}{a numeric vector}
    \item{\code{A2}}{a numeric vector}
    \item{\code{A3}}{a numeric vector}
    \item{\code{A4}}{a logical vector}
  }
}

\source{
Balding DJ, Buckleton J, Interpreting low template DNA profiles,
Forensic Science International: Genetics, 4: 1-10, 2009, doi: 10.1016/j.fsigen.2009.03.003.
http://www.zebfontaine.eclipse.co.uk/djb.htm
}
\references{
}
\examples{
data(Bates.DNA)
}
\keyword{datasets}