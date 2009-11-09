\encoding{UTF-8}
\name{simupopD}
\alias{simupopD}

\title{Simulate multi-population allele frequencies  for independent loci from a reference population,
 following a Dirichlet model}

\description{
Simulate multi-population allele frequencies for independent loci, from a given reference population, 
following a Dirichlet model.
Allele frequencies in the populations are generated as random deviates from
a Dirichlet distribution, whose parameters control the deviation of allele frequencies from the values
in the reference population. 
}

\usage{
simupopD(npop = 1, nloc = 1, na = 2, globalfreq = NULL, which.loc = NULL,
alpha1, alpha2 = 1)
}

\arguments{
  \item{npop}{ the number of populations}
  \item{nloc}{ the number of loci}
  \item{na}{ an integer vector giving the numbers of alleles per locus}
  \item{globalfreq}{ matrix of allele frequencies in the reference population. Data must be given in the format of the Journal of Forensic 
  Sciences for genetic data. Default corresponds to allele frequencies generated form a Dirichlet distribution 
  with parameter \code{alpha2} for all allele frequencies. }
  \item{which.loc}{ which loci to simulate from the \code{globalfreq} matrix, default considers all loci }
  \item{alpha1}{ a positive float  vector of length \code{npop} giving  the variance parameter  of the Dirichlet 
  distribution used
  to generate allele frequencies in the \code{npop} independent populations }
  \item{alpha2}{a positive float giving the parameter to be used to in the Dirichlet distribution to  generate 
  allele frequencies for the reference population }
}
\details{
 In the reference population, allele frequencies for independent loci  are simulated using a Dirichlet distribution with 
 parameter \code{alpha2}.\cr
 At  a given locus  L with n alleles,  the  allele frequencies are modeled as a vector of random 
 variables p=(p1, ..., pn) following a Dirichlet distribution with a  parameter vector of length n, 
 where each component is equal to alpha2,  p1+...+pn=1  and alpha2 > 0.\cr  
 Note that a  more sophisticated generation of global allele frequencies is possible using  the \code{\link{simufreqD}} function.
 Similarly, allele frequencies in the independent populations are simulated using a Dirichlet Distribution. 
 For example, for the first population to simulate, at a given locus  L with n alleles, 
 the  allele frequencies are modeled as a vector 
 of random  variables p=(p1, ..., pn) following a Dirichlet distribution with a  parameter vector of length n:\cr
 (p1(1-a1)/alpha1[1], ..., pn(1-alpha1[1])/alpha1[1]),  where p1+...+pn=1  and alpha1[1] > 0.\cr
alpha1[1] is the variance parameter for population 1 and is equivalent to Wright's Fst. The closest this parameter is to one,
the more the population allele frequencies are different from the values  of the reference population. 

}
\value{
  The result is stored in a list with two elements :
  \item{globfreq }{ a \code{tabfreq} object giving the allele frequencies of the chosen reference population,
  with the  chosen loci. }
  \item{popfreq}{a \code{tabfreq} object giving the allele frequencies of the simulated populations.}

}
\references{
Nicholson G, Smith AV, Jonsson F, Gustafsson O, Stefansson K, Donnelly P.
 Assessing population differentiation and isolation from single-nucleotide polymorphism data. 
 \emph{J Roy Stat Soc B} 2002;64:695–715 \cr
 
Marchini J, Cardon LR. Discussion on the meeting on "Statistical modelling and analysis of genetic data" 
 \emph{J Roy Stat Soc B}, 2002;64:740-741\cr
 
Wright S. The genetical structure of populations. \emph{Ann Eugen} 1951;15:323-354


 }

\author{Hinda Haned \email{haned@biomserv.univ-lyon1.fr}}
\note{ The code used here  for the generation of random Dirichlet deviates was previously implemented in the gtools library.}
\seealso{ \code{\link{simufreqD}}}
 
 
 	
 
\examples{
# simulate allele frequencies for two populations
data(Tu)
simupopD(npop=2,globalfreq=Tu, which.loc=c("FGA","TH01","TPOX"),
alpha1=c(0.2,0.3),alpha2=1)
  
}
\keyword{datagen}
