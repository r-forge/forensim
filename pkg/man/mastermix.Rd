\name{mastermix}
\Rdversion{1.1}
\alias{mastermix}
\title{
A Tcl/Tk graphical user interface for simple DNA mixtures resolution using allele peak heights/ or
 areas information
}
\description{
The \code{mastermix} function launches a Tcl/Tk graphical user interface dedicated 
to the resolution of two-person DNA mixtures using allele peak heights/ or areas information. 
\code{mastermix} is the implementation of a method developed by Gill et al (see the references section),
and previously programmed into an Excel macro by Dr. Peter Gill. 
}
\usage{
mastermix()
}



\details{
\code{mastermix} is  a Tcl/Tk graphical user interface implementing a  method developed by  Gill et al 
(1998) for simple mixtures resolution, using allele peak heights or areas information.  

This method searches through simulation the most likely combination(s) of the contributors' genotypes.
Having previously obtained an estimation for the mixture proportion,
it is possible to reduce the number of possible genotype combinations by keeping only 
those supported by the observed data. This is achieved by computing the sum of square differences between the expected 
allelic ratio and the observed allelic ratio, for all possible mixture combinations. 
The likelihood of peak heights (or areas), conditional on the combination of genotypes, is 
high if the residuals are low. 
Genotype combinations are thus selected according to the peak heights with the highest (conditioned) likelihoods. \cr

\code{mastermix} offers a graphical representation of the simulation for three models: 
\item The two allele model: at a given locus, two alleles are observed in the DNA stain
\item The three allele model: at a given locus, three alleles are observed in the DNA stain
\item The four allele model: at a given locus, four alleles are observed in the DNA stain

%Clicking on each button launches a simulation  dialog window for the corresponding model.
A left-click on each button launches a simulation dialog window for the corresponding
model, while a right-click opens the corresponding help page. 
}

\note{Each implemented model can either  be launched  using the  \code{mastermix} interface, or the 
\code{\link{A2.simu}}, \code{\link{A3.simu}} and  \code{\link{A4.simu}} functions, depending on the considered model.
}

\references{
Gill P, Sparkes P, Pinchin R, Clayton, Whitaker J, Buckleton J. Interpreting simple STR mixtures using allele peak
areas. \emph{Forensic Sci Int} 1998;91:41-5.
}
\author{ Hinda Haned  \email{haned@biomserv.univ-lyon1.fr} }
\seealso{\code{\link{A2.simu}}, \code{\link{A3.simu}} and {\code{\link{A4.simu}}}}

 

\examples{
mastermix()
}


\keyword{htest}