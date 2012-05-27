\name{PV}
\alias{PV}
\title{
Predictive value of the maximum likelihood estimator of the number of contributors to a DNA mixture
}
\description{
The \code{PV} function implements the predictive value of the maximum likelihood estimator of the number of contributors to a DNA
mixture
}
\usage{
PV(mat, prior)
}

\arguments{
  \item{mat}{matrix giving the estimates of the conditional probabilities that the maximum likelihood estimator
   classifies a given stain as a mixture of i contributors given that there are k contributor(s) to the stain. 
   Estimates i must be given in columns for each possible value of the number of contributors given in rows.  

}
  \item{prior}{ numeric vector giving the prior probabilities of encountering a mixture of i contributors. \code{prior} must be of 
  length the number of rows in \code{mat}. 
  }

}


\value{Vector of the predictive values}
\references{
Haned H., Pene L., Sauvage F., Pontier D., The predictive value of the maximum likelihood estimator of the number of contributors
to a DNA mixture, submitted, 2010.
}

\author{Hinda Haned  <h.haned@nfi.minvenj.nl>}

\seealso{
maximum likelihood estimator \code{\link{likestim}}
}


\examples{
# the following examples reproduce some of the calculations appearing
# in the article cited above, for illustrative purpose, the maximum 
#number of contributors is set here to 5 
#matcondi: Table 2 in Haned et al. (2010)
matcondi<-matrix(c(1,rep(0,4),0,0.998,0.005,0,0,0,0.002,0.937,0.067,0,0,0,0.058,
0.805,0.131,rep(0,3),0.127,0.662,rep(0,3),0.001,0.207),ncol=6)
#prior defined by a forensic expert (Table 3 in Haned et al., 2010)
prior1<-c(0.45,0.04,0.30,0.15,0.06)
#uniform prior, for each mixture type, the probability of occurrence is 1/5,
#5 being the threshold for the number of contributors
prior2<-c(rep(1/5,5))
#predictive values for prior1
PV(matcondi,prior1)
#for prior2
PV(matcondi, prior2)
}

\keyword{htest}
