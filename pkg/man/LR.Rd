\name{LR}
\alias{LR}

\title{ Likelihood ratio for DNA evidence interpretation}
\description{
The \code{LR} function calculates the likelihood ratio for a DNA evidence, when two 
competing hypotheses Hd and Hp, respectively the defence and the prosecution hypotheses,
are weighted about the origin of the DNA evidence. The evidence can either be a simple
or a mixed stain.  

}
\usage{
LR(stain, freq, xp=0, xd=0, Tp=NULL, Vp=NULL, Td=NULL, Vd=NULL, theta=0)
}
\arguments{
  \item{stain}{
a vector giving the set of (distinct) alleles present in the DNA stain
}
  \item{freq}{
vector of the corresponding allele frequencies in the global population
}
  \item{xp}{
  the number of unknown contributors to the stain under the prosecution hypothesis Hp. Default is 0.
}
  \item{xd}{
 the number of unknown contributors to the stain under the defence hypothesis Hd. Default is 0.
}
\item{Tp}{
a vector of strings where each string contains two alleles separated by '/', 
corresponding to one known contributor under the prosecution hypothesis Hp. 
The length of the vector equals the number of known contributors. Default is NULL.
}
\item{Vp}{
a vector of strings where each string contains two alleles separated by '/', 
corresponding to one known non-contributor under the prosecution hypothesis Hp. The length of 
the vector equals the number of known non-contributors.  Default is NULL.
}
\item{Td}{
a vector of strings where each string contains two alleles separated by '/', 
corresponding to one known contributor under the defence hypothesis Hd. The length of the
vector equals the number of known contributors.  Default is NULL.
}
  \item{Vd}{
a vector of strings where each string contains two alleles separated by '/', 
corresponding to one known non-contributor under the defence hypothesis Hd. The length of the
 vector equals the number of known non-contributors.  Default is NULL.
}

\item{theta}{ a float  in [0,1[. \code{theta} is equivalent to Wright's Fst.
In case of population subdivision, it allows a correction of the allele frequencies in the
subpopulation of interest}
}




\details{
\code{LR} is the implementation of the general formula of Curran et al (1999) for the
evaluation of forensic DNA mixtures through likelihood ratios. The likelihood ratio is computed 
as a ratio of two probabilities of the DNA evidence, \eqn{E}, conditional on the evaluated hypotheses:
\deqn{LR=\frac{P(E|H_p)}{P(E|H_d)},}{P(E|H_p) / P(E|H_d),}
where \eqn{H_p} denotes the prosecution hypothesis and \eqn{H_d} the defence hypothesis. 

In case of population subdivision, contributors to the DNA stain are considered to come from the same 
subpopulation.
Allele dependencies within subpopulations are accounted for through Wright's Fst coefficient, denoted here
\eqn{\theta}.
}

\references{
Curran JM, Triggs CM, Buckleton J, Weir BS. Interpreting DNA Mixtures in Structured Populations. 
\emph{ J Forensic Sci} 1999;44(5): 987-995

}
\author{Hinda Haned \email{h.haned@nfi.minjus.nl}}

\note{
Please note that the \code{LR} function is based on functions initially implemented in
the forensic package by Miriam Marusiakova
 \url{http://cran.r-project.org/web/packages/forensic/}}




\seealso{
 the exclusion probability \code{\link{PE}}.}

\examples{
# A rape case in Hong Kong (Hu and Fung, Int J Legal Med 2003)
# The stain shows alleles 14, 15, 17 and 18 at locus D3S1358.
stain =c(14,15,17,18)
# suspect's profile: "14/17"
suspect<-"14/17"
# victim's profile: "15/18"
victim<-"15/18"
# corresponding allele frequencies 
freq<-c(0.033,0.331,0.239,0.056)

# Prosecution hypothesis: Contributors were the victim and the suspect 
# defence hypothesis: Contributors were the victim and 1 unknown contributor
# Likelihood ratios for DNA evidence for different alternatives:
LR(stain,freq,xp=0,Tp=c(victim,suspect),Vp=NULL,Td=victim,Vd=suspect,xd=1)

}

\keyword{htest}
