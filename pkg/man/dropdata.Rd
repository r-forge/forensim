\encoding{UTF-8}
\name{dropdata}
\alias{dropdata}
\docType{data}
\title{Dropout example data}
\description{
\code{dropdata} gives is an extract of a series of experiments used to determine the probability of dropout
}

\usage{data(dropdata)}
\format{
A list of two components: 'tabcsv' and 'genot'.
\code{tabcsv} is an extract of the validation table of a two-person mixture (Genemapper format)  and
\code{genot} is the matrix of genotypes  of the  individuals contributing to the mixture.
}

\details{
The mixture is characterized using the Applied Biosystems AmpFlSTR Identifiler^{TM} kit.
}

\source{
Data communicated by Elodie Suzanne and Laurent P\'ene, Laboratoire de Police Scientifique,
Ecully, France.
}
\examples{
data(dropdata)
names(dropdata)
dropdata$tabcsv
dropdata$genot

}
\keyword{datasets}
 