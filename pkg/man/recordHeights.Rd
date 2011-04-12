\name{recordHeights}
\alias{recordHeights}
\title{
Records the peak heights of the alleles present in the analyzed stains
}
\description{
The \code{recordHeights} function records the peak heights of the alleles present in the analyzed stains.
The function aims to facilitate the 
manipulation of experimental data used for the estimation of the probability of allelic dropout
(cf. the references sections).
}
\usage{
recordHeights(x,y=NULL,geno,tabcsv,byloc=FALSE)
}

\arguments{
  \item{x}{
numeric label of the contributing individual, if the stain is a mixture, \code{x} should give the
label of the first individual contributing to the mixture
}
  \item{y}{
numeric label of the second contributing individual, default is NULL. If the stain is a mixture, \code{y}
should give the label of the second individual contribution to the mixture. 
This argument is skipped if the stain is not a mixture (default case: \code{y} set to NULL).
}
  \item{geno}{
a matrix giving the genotypes of the individuals contributing to the analyzed data for each locus. 
An individual genotype is given in rows for each locus in column. A homozygous carrying allele 9 is coded
'9/9', a heterozygous carrying alleles 8 and 9 is coded '8/9'. Individual labels are coded using integers that
are simply the order of introduction in the data frame.

}
  \item{tabcsv}{
a matrix giving the validation table of the analysed DNA stain. \code{tabcsv} must have a ``genemapper'' validation
table structure, namely, information about the present alleles and the corresponding peak heights must be given. 
}
  \item{byloc}{
logical indicating whether data should be displayed per locus (TRUE) or overall loci (FALSE, default)
}
}


\value{
A list of length the number of analyzed loci, each component of the list is a matrix with the following information:
- The names of expected alleles
- The expected allele counts for the first contributor (when date is a mixture)
- The expected allele counts for the second contributor (when date is a mixture)
- The observed alleles
- The observed peak heights
- The dropout variable D, takes 1 if the allele has dropped out, 0 otherwise


}
\references{
Gill P, Puch-Solis R, Curran J. The low-template-DNA (stochastic) threshold-Its determination relative to risk 
analysis for national DNA databases. \emph{Forensic Science International: Genetics}, 2009, 3, 104-111
}
\author{Hinda Haned  <h.haned@nfi.minjus.nl>}
\note{

\code{recordHeights} is an auxiliary function of the \code{\link{tabSPH}} function that implements the 
methodology proposed by Gill et al. to estimate the probability of allelic dropout using experimental DNA mixtures.
}

\seealso{\code{\link{recordDrop}} for an alternative method, \code{\link{tabSPH}}
}
\examples{
#load the exemple data
data(dropdata)
tabcsv<-dropdata$tabcsv
genot<-dropdata$genot
#individuals' labels are 1 and 2
#peak heights of heterozygote genotypes with non shared alleles
recordHeights(1,2,geno=genot,tabcsv=tabcsv)
}
\keyword{manip}
\keyword{datagen}