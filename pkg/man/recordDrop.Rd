\name{recordDrop}
\alias{recordDrop}
\title{
Records the allelic dropout events matched with individual DNA proxies }
\description{
The \code{recordDrop} function records the dropout events from experimental data. The function aims to facilitate the 
manipulation of experimental data used for the estimation of the probability of allelic dropout (cf. the references sections).
}
\usage{
recordDrop(x, y, geno, tabcsv,s=40)
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
\item{s}{
numeric giving the detection threshold for alleles in Relative fluorescence units (RFU) , default is set to 40 RFUS. An observed
allele with a peak height smaller (<) than 40 RFUS is considered as dropped-out. 
}
}

\value{
A list of length the number of analyzed loci, each component of the list is a matrix with the following information:
- The names of expected alleles\cr
- The expected allele counts for the first contributor (when date is a mixture)\cr
- The expected allele counts for the second contributor (when date is a mixture)\cr
- The observed alleles\cr
- The observed peak heights\cr
- The dropout variable D, takes 1 if the allele has dropped out, 0 otherwise\cr


}
\references{
Tvedebrink T, Eriksen PS, Mogensen HS, Morling N. Estimating the probability of allelic drop-out of STR alleles
in forensic genetics. \emph{Forensic Science International: Genetics}, 2009, 3(4), 222-226.
}
\author{Hinda Haned  <haned@biomserv.univ-lyon1.fr>}
\note{

\code{recordDrop} is an auxiliary function of the \code{\link{tabDNAproxy}} function that implements the 
methodology proposed by Tvedebrink et al. to estimate the probability of allelic dropout using experimental DNA mixtures.
}

\seealso{\code{\link{DNAproxy}}, \code{\link{tabDNAproxy}}
}
\examples{
#load the exemple data
data(dropdata)
tabcsv<-dropdata$tabcsv
genot<-dropdata$genot
#individuals' labels are 1 and 2
#record the dropout the surviving peak heights for heterozygotes whith non shared alleles
recordDrop(1,2,geno=genot,tabcsv=tabcsv,s=40)

}
\keyword{manip}
\keyword{datagen}