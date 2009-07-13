\name{A2.simu}
\Rdversion{1.1}
\alias{A2.simu}
\title{

A Tcl/Tk graphical user interface for simple DNA mixtures resolution using allele peak heights or areas information 
when two alleles are observed at a given locus
}
\description{
The \code{A2.simu} function launches a Tcl/Tk graphical interface with functionalities devoted to two-person DNA mixtures 
resolution, when  two alleles are observed at a given locus. 

}
\usage{
A2.simu()
}



\details{
When two alleles are observed at a given locus  in the DNA stain, seven genotype combinations
are possible for the two contributors: (AA,AB), (AB,AB), (AA,BB), (AB,AA), (BB,AA), (AB,BB) and (BB,AB), where A and B are 
the two observed alleles (in ascending order of molecular weight).
 Having previously obtained an estimation for the mixture proportion,
it is possible to reduce the number of possible genotype combinations by keeping those only supported by the 
observed data. This is achieved by computing the sum of square differences between the expected 
allelic ratio and the observed allelic ratio, for all possible mixture combinations. 
The likelihood of peak heights (or areas), given the combination of genotypes, is 
high if the residuals are low. 
Genotype combinations are thus selected according to the peak heights with the highest likelihoods. \cr

The \code{A2.simu()} function launches a dialog window with three buttons:
\item \code{Plot simulations}: plot of the residuals of each possible genotype combination for varying values of the mixture proportion across the interval [0.1, 0.9].
The observed mixture proportion is also reported on the plot. 
\item \code{Simulation details}: a matrix containing the simulation results. Simulation details and genotype combinations
with the lowest residuals can be saved as a text file by clicking the 
``Save" button.  It is also possible to choose specific paths and names for the save files
\item \code{Genotypes filter}: a matrix giving the mixture proportion conditional on the genotype combination. This conditional mixture proportion helps filter the most
plausible genotypes among the seven possible combinations. The matrix can be saved as a text file by clicking the ``Save" button. 
It is also possible to choose a specific path and a name for the save file.\cr
}


\note{
\item Please note that some buttons won't be working the same according to your platform, this is because
the tcl/tk package (more specifically, the tktable widget) is not up to date under Linux and MacOS systems. 
Consequently, under these systems, the ``Simulations details" and ``Genotypes filter" buttons yield an interface
to save the results, and the corresponding tables can only be displayed under the Windows system.


\item Please see forensim tutorial for the computational details at \url{http://forensim.r-forge.r-project.org/misc/forensim-tutorial.pdf}
}
\references{
Gill P, Sparkes P, Pinchin R, Clayton, Whitaker J, Buckleton J. Interpreting simple STR mixtures using allele peak
areas. \emph{Forensic Sci Int} 1998;91:41-53.
}
\author{ Hinda Haned  \email{haned@biomserv.univ-lyon1.fr} }



 \seealso{
\code{\link{A3.simu}}: the three-allele model, and \code{\link{A4.simu}}: the four-allele model}
 
\examples{
A2.simu()
}


\keyword{htest}