\encoding{UTF-8}
\name{simPCR2}
\alias{simPCR2}
\title{Polymorphism chain reaction simulation model}
\description{\code{simPCR2} implements a simulation model for the polymorphism chain reaction (Gill et al., 2005).
Giving several input parameters, \code{simPCR2}  outputs the number of amplified DNA molecules and their corresponding peak heights
(in RFUs).
 }

 
\usage{
simPCR2(ncells,probEx,probAlq, probPCR, cyc = 28, Tdrop = 2 * 10^7,
probSperm = 0.5, dip = TRUE,KH=55)
}
\arguments{
\item{ncells}{initial number of cells}
\item{probEx}{ probability that a DNA molecule is extracted (probability of surviving the extraction process) }
\item{probAlq}{probability that a DNA molecule is selected for PCR amplification}
\item{probPCR}{probability that a DNA molecule is amplified during a given round of PCR}
\item{cyc}{number of PCR cycles, default is 28 cycles}
\item{Tdrop}{threshold of detection: number of molecules (in the total PCR reaction mixture)
 that is needed to generate a signal, default is set to 2*10^7 molecules}
\item{probSperm}{probability of observing alleles of type A in the initial sample of haploid cells (e.g. sperm cells). Probability
of observing allele B is given by 1-\code{probSperm}}
\item{dip}{ logical indicating the cell ploidy, default is diploid cells (TRUE), FALSE is for haploid cells}
\item{KH}{positive constant used to scale the peak heights obtained from the number of amplified molecules
(see reference section)}
}


\details{
A threshold of \code{Tdrop} (must be a multiple of 10^7) is needed to generate a signal, then, a log-linear relationship is used to
determine the intensity of the signal with respect to the number of successfully amplified DNA molecules. Dropout events occur
 whenever less than \code{Tdorp} molecules are generated.
}

\value{
A matrix with the following components:
\item{HeightA}{Peak height of allele A}
\item{DropA}{Dropout variable for allele A}
\item{HeightB}{Peak height of allele B}
\item{DropB}{Dropout variable for allele B}
}

\references{

Jeffreys AJ, Wilson V, Neumann R and Keyte J.  Amplification of human minisatellites by the polymerase chain reaction: towards
DNA fingerprinting of single cells. \emph{Nucleic Acids Res} 1988;16: 10953_10971.\cr

Gill P, Curran J and Elliot K. A graphical simulation model of the entire DNA process associated with the analysis of short 
tandem repeat loci.\emph{ Nucleic Acids Research} 2005, 33(2): 632-643.
}

\seealso{ \code{\link{simPCR2TK}}}  

\author{ Hinda Haned  \email{haned@biomserv.univ-lyon1.fr} }



\examples{

#simulation of a 28 cycles PCR, with the initial stain containing 5 cells 
simPCR2(ncells=5,probEx=0.6,probAlq=0.30,probPCR=0.8,cyc=28, Tdrop=2*10^7,dip=TRUE,KH=55)

}
%\keyword{manip}
\keyword{datagen}

