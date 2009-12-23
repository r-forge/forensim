\name{simMixSNP}
\alias{simMixSNP}
\title{
Simulates SNP mixtures
}
\description{
Simulates SNP mixtures and outputs optionally file
suitable for wrapdataL function
for estimation of number of contributors
}
\usage{
simMixSNP(nSNP , p , ncont, writeFile, outfile , id )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nSNP}{
Integer number of SNPs>1
}
  \item{p}{
Minor allele frequency
}
  \item{ncont}{
Number of contributors >= 1
}
  \item{writeFile}{
If TRUE, output written to file
}
  \item{outfile}{
Name of output file
}
  \item{id}{
Column one of output file identifying run
}
}
\value{
Returns a data frame with columns Id, marker, allele, frequency and height (=1 for now)
}
\author{
Thore Egeland <Thore.Egeland@medisin.uio.no>
}

\examples{
simMixSNP()
}

