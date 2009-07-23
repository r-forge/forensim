\name{Pevid2}
\alias{Pevid2}
\title{ Conditional profile probabilities}

\description{
Calculates the probability of observing a set of DNA profiles conditional on a given hypothesis
specifying who were the contributors to the observed profiles.
All the  individuals involved in the analyzed case are assumed to come from the same
subpopulation with a given coancestry coefficient. 
}
\usage{
Pevid2(stain, freq, x, T = NULL, V = NULL, theta = 0)
}
\arguments{
  \item{stain}{ vector of distinct alleles (from one specific locus) found in 
   the crime sample.}
  
  \item{freq}{vector of the corresponding allele frequencies in the global population}
  
  \item{x}{ the number of unknown contributors to the mixture }
  
  \item{T}{ object of class genotype (package \pkg{genetics}), or a vector of strings 
    where each string contains two alleles separated by '/', corresponding to one known 
    contributor. The length of the vector equals the number of known contributors. 
    Default is NULL.
				 }
  \item{V}{ object of class genotype (package \pkg{genetics}), or a vector of strings 
    where each string contains two alleles separated by '/', corresponding to one known 
    non-contributor. The length of the vector equals the number of known non-contributors. 
    Default is NULL.    
 }
  \item{theta}{ a float  in [0,1[. \code{theta} is equivalent to Wright's Fst. In case of
  population
  subdivision, it allows a correction of the allele frequencies in the subpopulation of 
  interest}
 }


\references{ 
    Curran JM, Triggs CM, Buckleton J, Weir BS. Interpreting DNA Mixtures in Structured Populations. 
\emph{ J Forensic Sci} 1999;44(5): 987-995
     
 }
\author{Hinda Haned \email{haned@biomserv.univ-lyon1.fr}}

\note{
Please note that the \code{Pevid2} function is an improved version of the \code{Pevid.gen}
function from the forensic package by Miriam Marusiakova (which explains the 2 in the function
name). \code{Pevid2} calls external functions in C code. 
}


\seealso{\code{\link{LR}}, \code{\link{RMP}}} 
\examples{
# A rape case in Hong Kong (Hu and Fung, Int J Legal Med 2003)
# The stain shows alleles 14, 15, 17 and 18 at locus D3S1358.
stain=c(14,15,17,18)
# suspect's profile: "14/17"
suspect<-"14/17"
# victim's profile: "15/18"
victim<-"15/18"
# corresponding allele frequencies 
freq<-c(0.033,0.331,0.239,0.056)

# Prosecution proposition: Contributors were the victim and the suspect 
# defence proposition: Contributors were the victim and 1 unknown contributor
# from the same subpopulationas the victim
# Evaluation of the defence proposition, in case of independence between alleles
Pevid2(stain, freq, x=1, T = victim)

# note that if theta=0, the suspect's profile plays no role in the calculation
#and the same result is obtained
Pevid2(stain, freq, x=1, T = victim, V = suspect)
# In case of allele dependencies, measured by theta=0.03
Pevid2(stain, freq, x=1, T = victim, V = suspect, theta = 0.03)

}

\keyword{ htest }
