\name{tabDNAproxy}
\alias{tabDNAproxy}

\title{
Builds a list of tables that record the dropout events matched with the appropriate DNA proxies
}
\description{
The \code{tabDNAproxy} function builds a list of tables that record the dropout events matched with the appropriate ``DNAproxies'', these are the approximations
of the amount of DNA contributed by the individuals in the analyzed DNA stains. Each table is specific to a locus. This function builds the data frames on which
the logistic model, proposed by Tvedebrink et al (cf. references section), can be performed.
}
\usage{
tabDNAproxy(x, y = NULL, geno, tabcsv)
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

}

\value{
A list of length the number of analyzed loci, each component of the list is a matrix with the following information:
\item{Dloc}{the (per locus) dropout variable D, takes 1 if the allele has dropped out, 0 otherwise}
\item{Hestim }{mean peak heights derived from the DNA proxies, see the references section for further details }
}
\references{
Tvedebrink T, Eriksen PS, Mogensen HS, Morling N. Estimating the probability of allelic drop-out of STR alleles
in forensic genetics. \emph{Forensic Science International: Genetics}, 2009, 3(4), 222-226.
}
\author{Hinda Haned  <haned@biomserv.univ-lyon1.fr>}



\seealso{
\code{\link{recordDrop}}, \code{\link{DNAproxy}}
}
\examples{
#load the exemple data
data(dropdata)

tabcsv<-dropdata$tabcsv
genot<-dropdata$genot
#individuals' labels are 1 and 2
#lets record the dropout events and the corresponding DNA proxies
tabDNAproxy(1,2,geno=genot,tabcsv=tabcsv)

}
\keyword{manip}
\keyword{datagen}

