\name{A3.simu}
\alias{A3.simu}
\title{

A Tcl/Tk graphical user interface for simple DNA mixtures resolution using allele peak heights or areas when three alleles are observed at a given locus
}
\description{
The \code{A3.simu} function launches a Tcl/Tk graphical interface with functionalities devoted to two-person
 DNA mixtures resolution, when three alleles are observed at a given locus. 

}
\usage{
A3.simu()
}



\details{
When  three alleles are observed at a given locus  in the DNA stain, twelve genotype combinations
are possible for the two contributors: (AA,BC), (BB,AC), (CC,AB), (AB,AC), (BC,AC), (AB,BC), (BC,AA), (AC,BB), (AB,CC), (AC,AB), (AC,BC) and (BC,AB) where A, B and C are the three
observed alleles (in ascending order of molecular weights). 
Having previously obtained an estimation for the mixture proportion,
it is possible to reduce the number of possible genotype combinations by keeping those only supported by the 
observed data. This is achieved by computing the sum of square differences between the expected 
allelic ratio and the observed allelic ratio, for all possible mixture combinations. 
The likelihood of peak heights (or areas), given the combination of genotypes, is 
high if the residuals are low. 
Genotype combinations are thus selected according to the peak heights with the highest  likelihoods. \cr

The \code{A3.simu()} function launches a dialog window with three  buttons:\cr
-\code{Plot simulations}: plot of the residuals of each possible genotype combination for varying values of the mixture proportion across the interval [0.1, 0.9].
The observed mixture proportion is also reported on the plot.\cr
-\code{Simulation details}: a matrix containing the simulation results. Simulation details and genotype combinations
with the lowest residuals can be saved as a text file by clicking the 
``Save" button.  It is also possible to choose specific paths and names for the save files.\cr
-\code{Genotypes filter}: a matrix giving the mixture proportion conditional on the genotype combination. This conditional mixture proportion helps filter the most
plausible genotypes among the twelve possible combinations. The matrix can be saved as a text file by clicking the ``Save" button. 
It is also possible to choose a specific path and a name for the save file.\cr
}

\note{

-Linux users may have to download the \code{libtktable} package to their system before using the \code{A3.simu} function. 
This is due to the \code{Tktable} widget,  used in \code{forensim}, which is not (always)
downloaded with the Tcl/Tk package.\cr
-For the computational details, please see forensim tutorial at \url{http://forensim.r-forge.r-project.org/misc/forensim-tutorial.pdf}.
}

\references{
Gill P, Sparkes P, Pinchin R, Clayton, Whitaker J, Buckleton J. Interpreting simple STR mixtures using allele peak
areas. \emph{Forensic Sci Int} 1998;91:41-53.
}
\author{ Hinda Haned  \email{haned@biomserv.univ-lyon1.fr} }



 \seealso{
\code{\link{A2.simu}}: the two-allele model, and \code{\link{A4.simu}}: the four-allele model}
 
\examples{
A3.simu()
}

\keyword{htest}