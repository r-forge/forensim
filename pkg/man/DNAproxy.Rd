\name{DNAproxy}
\alias{DNAproxy}
\title{
Approximation of the amount of DNA contributed by a person based on the observed peak heights of the alleles present in the analyzed sample
}
\description{
\code{DNAproxy} gives an estimation of the amount of DNA contributed by a person to a DNA stain based on the observed peak 
heights of the present alleles. The estimation is performed using data across all available loci, data can either consist of
single-contributor or mixed DNA stains. The computation of the DNA proxies from experimental data are
described by Tvedebrink et al. (cf. the references sections).}
\usage{
DNAproxy(tab, x)
}
\arguments{
  \item{tab}{
a table produced by the \code{\link{recordDrop}} function, giving the allelic dropouts observations and the corresponding
allelic peak heights
}
  \item{x}{
a character giving the label of the individual for whom the DNA proxy must be specified, this argument is to be specified only
when data in \code{tab} is made of mixtures. In case data is consist of single-contributor stains, the argument must be left empty.
\code{x} must match the name given in the \code{tab} table.
}
}

\references{
Tvedebrink T, Eriksen PS, Mogensen HS, Morling N. Estimating the probability of allelic drop-out of STR alleles
in forensic genetics. \emph{Forensic Science International: Genetics}, 2009, 3(4), 222-226.
}
\author{Hinda Haned  <h.haned@nfi.minjus.nl>}

\note{

\code{DNAproxy} is an auxiliary function of the \code{\link{tabDNAproxy}} function that implements the 
methodology proposed by Tvedebrink et al. to estimate the probability of allelic dropout using experimental DNA mixtures.
}
\seealso{
\code{\link{recordDrop}}, \code{\link{tabDNAproxy}}
}
\examples{
#load the exemple data
data(dropdata)
tabcsv<-dropdata$tabcsv
genot<-dropdata$genot
#individuals' labels are 1 and 2
#DNA proxy for individual one, when data is composed of a 2-person mixture
DNAproxy(recordDrop(1,2,geno=genot,tabcsv=tabcsv),"c1")

}

\keyword{manip}
\keyword{datagen}

