\name{tabSPH}
\alias{tabSPH}

\title{
Builds a matrix of the dropout variable and the corresponding surviving peak heights
}
\description{
The \code{tabSPH} function builds a matrix of the dropout variable and the corresponding surviving peak heights,
for each available locus or across all loci (default).  The constructed matrices have two columns: the dropout
variable and the surviving peak heights. The logistic model proposed to model the dropout probability from
experimental data (see the references section) can be performed directly on the data yielded by \code{tabSPH}.  
}
\usage{
tabSPH(x, y = NULL, geno, tabcsv, byloc = FALSE,s=40)
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

\item{s}{
numeric giving the detection threshold for alleles in Relative fluorescence units (RFU) , default is set to 40 RFUS. An observed
allele with a peak height smaller (<) than 40 RFUS is considered as dropped-out. 
}
}
\details{
Both mixed and unmixed samples can be used in \code{tabSPH}, setting  the \code{y} argument to NULL (default)
will produce results considering data for \code{x} only. In case of mixtures, note that only heterozygote
genotypes with no shared  alleles are considered.
}
\value{
If argument \code{byloc} is TRUE, \code{tabSPH} yields a list of length the number of available loci, each elements of the list
contain a matrix with two columns:
\item{D}{the dropout variable}
\item{H}{the surviving peak height}
If argument \code{byloc} is FALSE, \code{tabSPH} yields a single matrix with columns D and H described  above.
 }
\references{
Gill P, Puch-Solis R, Curran J. The low-template-DNA (stochastic) threshold-Its determination relative to risk 
analysis for national DNA databases. \emph{Forensic Science International: Genetics}, 2009, 3, 104-111

}
\author{Hinda Haned  <haned@biomserv.univ-lyon1.fr>}



\seealso{
\code{\link{tabDNAproxy}}
}
\examples{
#load the example data
data(dropdata)
tabcsv<-dropdata$tabcsv
genot<-dropdata$genot
#individuals' labels are 1 and 2
#recording dropout variable matched with the surviving peak heights
#for heterozygotes whith non shared alleles
tabSPH(1,2,geno=genot,tabcsv=tabcsv,s=0)

}

\keyword{manip}
\keyword{datagen}