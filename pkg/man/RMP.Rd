\name{RMP}
\alias{RMP}

\title{
The Random Match Probability of DNA evidence (RMP)
}
\description{
\code{RMP} computes the random match probability of DNA evidence given in a matrix (or data frame) or
in a text file. Several situations are handled: the suspect and an unknown offender are unrelated,
or are members of the same subpopulation with a given coancestry coefficient theta, or are close
relatives. For the latter case, the relationship is described by the kinship coefficients.
 
}
\usage{
RMP(suspect=NULL, filename=NULL, freq, k=c(1,0,0), theta=0,refpop=NULL)
}

\arguments{
\item{suspect}{
a matrix or a data frame of dimension L x 2, L being the number of loci involved in the DNA
 evidence.
The first column gives the loci names, and the second column gives the suspect's genotype 
at each locus. A genotype is coded as a character where each string contains two alleles
 separated by '/'. 
The DNA evidence can also be given in a text file, see argument \code{filename}.
}


  \item{filename}{
the file name from which the input data should be read. Data mut be a matrix of dimension L x 2, L being the number of loci involved in the DNA
 evidence.The first column gives the loci names, and the second column gives the suspect's genotype 
at each locus. A genotype is coded as a character where each string contains two alleles
 separated by '/'. 
}
  \item{freq}{
a \code{tabfreq} object giving the allele frequencies 
}
  \item{k}{
vector of kinship coefficients (k0, k1, k2), where ki is the probability that two
people (the suspect and an unknown offender) will share i alleles identical by
descent, i = 0, 1, 2.
}
\item{theta}{ a float  in [0,1[. \code{theta} is equivalent to Wright's Fst.
In case of population subdivision, it allows a correction of the allele frequencies in the
subpopulation of interest}

 \item{refpop}{
the reference population in \code{freq} from which to extract the allele frequencies fro the RMP calculation.
This argument is obligatory only if \code{freq} contains allele frequencies from several populations 
}
}


\details{
The match probability is derived from Balding and Nichols (1994) and is computed as:
\deqn{k_2 + k_1 Z_1 + k_0 Z_2} where
\eqn{k_0, k_1, k_2} are the kinship coefficients,\cr
\eqn{Z_1} is the match probability when the suspect an the unknown offender share one allele identical-by-descent.\cr
\eqn{Z_2} is the match probability in the unrelated case, when the suspect an the unknown offender share 0 allele 
identical-by-descent.


In the homozygous case, with the allele frequency \eqn{p_i}:
\deqn{
Z_1=\frac{2\theta + (1-\theta)p_i}{1+\theta}}{%
Z_1 = (2\theta + (1-\theta) p_i)/ (1+\theta)
}
\deqn{
Z_2=\frac{\left[2\theta+(1-\theta)p_i\right] \left[3\theta+(1-\theta)p_i \right]}{(1+\theta)(1+2\theta)}}{%
Z_2 = [2theta + (1-theta)p_i][3theta + (1-theta)p_i] / [(1 + theta)(1 + 2theta)]}

In the heterozygous case, with allele frequencies \eqn{p_i} and \eqn{p_j}:

\deqn{Z_1=\frac{2\theta + (1-\theta)(p_i+p_j)}{2(1+\theta)}}{%
Z_1 = (2\theta + (1-\theta) (p_i+p_j))/ (2(1+\theta))
}
\deqn{Z_2=\frac{2\left[\theta+(1-\theta)p_i\right] \left[\theta+(1-\theta)p_j \right]}{(1+\theta)(1+2\theta)}}{%
Z_2 = 2[theta + (1-theta)p_i][theta + (1-theta)p_j] / [(1 + theta)(1 + 2theta)]}

 \eqn{\theta} is  Wright's Fst coefficient, usually called the coancestry coefficient in forensic studies.
Main effects of allele dependencies between loci in the suspect's subpopulation are taken into account 
though the coancestry coefficient, hence, the match probability at all loci is, to a close approximation,
the product of single-locus probabilities.
}
\value{
\code{RMP} returns a list with the following components:
\item{RMP.loc }{single-locus match probabilities}
\item{RMP}{multiloci match probability (product of single-locus match probabilities)}
}

\references{ 
Balding DJ, Nichols RA. DNA profile match probability calculation: How 
to allow for population stratification, relatedness, database selection and 
single bands. \emph{Forensic Sci I} 1994;64:125-140.
}

\author{Hinda Haned \email{haned@biomserv.univ-lyon1.fr}}


\seealso{ \code{\link{LR}} for the evaluation of DNA evidence through likelihood ratio}
\examples{
# random match probability 
# data input 

data <- matrix(c("CSF1PO","FGA","TH01","TPOX","VWA","D3S1358","D5S818",
"D7S820","D8S1179","D13S317","D16S539","D18S51","D21S11","D2S1338","D19S433",
"12/11","22/19","6/7","10/8","17/18","18/17","12/12","8/8","13/13","11/11",
"12/10","14/15","33.2/32.2","23/22","14/14"),nc=2)
colnames(data)<- c('locus','genotype')
#15-locus genotype
data
#allele frequencies are taken from the strusa data set

data(strusa)

RMP(suspect=data,freq=strusa,refpop="Cauc")

# using a preexisting file from the forensim package
RMP(filename=system.file("files/exprofile.txt", package = "forensim"),
freq=strusa,refpop="Cauc")
 
}
\keyword{ htest}
