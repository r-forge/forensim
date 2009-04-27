\encoding{UTF-8}
\name{strusa}
\alias{strusa}
\docType{data}
\title{  Allele frequencies for 15 autosomal short tandem repeats  core loci on U.S. Caucasian, African American, and Hispanic populations.}
\description{Allele frequencies for 15 autosomal short tandem repeats  loci on three American populations : Caucasians, African Americans and Hispanics.
Among the 15 loci, 13 belong to the core Combined DNA Index System (CODIS) loci used by the Federal Bureau of Investigation (USA),  in forensic DNA analysis, and two 
supplementary loci are more commonly used in Europe, see details.
 
 
}


\details{
CSF1PO, FGA, TH01, TPOX, vWA, D3S1358, D5S818, D7S820,
D8S1179, D13S317, D16S539, D18S51 and  D21S11, belong to the core CODIS loci used in the US, whereas D2S1338 and  D19S433 belong to the European core loci. 
}
\usage{data(strusa)}

\format{
    \code{strusa} is a tabfreq  object giving allele frequencies of 15 lociin three American populations. 
}

\references{
Butler JM, Reeder DJ.  http://www.cstl.nist.gov/strbase/index.htm, last visited: March 4th 2009 \cr
Butler JM, Schoske R, Vallone MP, Redman JW, Kline MC.  Allele frequencies for 15 autosomal STR loci on U.S. Caucasian, 
African American, and Hispanic populations. \emph{ J Forensic Sci} 2003;48(8):908-911.
}
\examples{
data(strusa)
strusa
#genotypes simulations from each population
geno<- simugeno(strusa,n=c(100,100,100))
geno
#3-person mixture simulation with the contributors from the three different populations
mix3<- simumix(geno,ncontri=c(1,1,1))
mix3

}
\keyword{datasets}
