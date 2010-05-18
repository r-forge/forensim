%\encoding{UTF-8}
\name{Tu}
\alias{Tu}
\docType{data}
\title{Allele frequencies of 15 autosomal short tandem repeats loci on Chinese Tu ethnic minority group}
\description{
Population genetic analysis of 15 STR loci of Chinese Tu ethnic minority group.
}

\usage{data(Tu)}
\format{

a data frame  presented in the format of the Journal of Forensic Sciences for genetic data: allele names are given in the first column, and frequencies for a given 
allele are read in rows for the different markers. When a given allele is not observed, value is coded NA 
(rather than "-" in the original format).
  

}

\details{
CSF1PO, FGA, TH01, TPOX, vWA, D3S1358, D5S818, D7S820,
D8S1179, D13S317, D16S539, D18S51 and  D21S11, belong to the core CODIS loci used in the US, whereas D2S1338 and  D19S433 belong to the European  core loci.
}

\references{
 Zhu B,  Yan J,  Shen C, Li T,  Li Y, Yu X,  Xiong X,  Muf H,  Huang Y, Deng Y. (2008).
 Population genetic analysis of 15 STR loci of Chinese Tu ethnic minority group. \emph{Forensic Sci Int}; 174: 255-258.
}
\examples{
data(Tu)
tabfreq(Tu)

}
\keyword{datasets}
 