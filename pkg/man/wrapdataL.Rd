\name{wrapdataL}
\alias{wrapdataL}
\title{
ML estimate of number of contributors for SNPs
}
\description{
Wrap up of dataL in forensim. Given file with columns:
"No, Marker, Allele, Frequency and  Height"
the log likelihood for requested number of contributors
is calculated. For now only "Frequency" column is used.
}
\usage{
wrapdataL(fil , plotte , nInMixture , tit )
}
\arguments{
  \item{fil}{
Input file
}
  \item{plotte}{
If T, plot
}
  \item{nInMixture}{
Alternatives for number of contributors, say 1:5
}
  \item{tit}{
Title to be used in plot
}
}
\value{
Plot (optional) and  log likelihoods
}
\author{
Thore Egeland \email{Thore.Egeland@medisin.uio.no}
}
\examples{
aa<-simMixSNP(nSNP=5,writeFile=TRUE,outfile="sim.txt",ncont=3) #Simulates data
res<-wrapdataL(fil="sim.txt") # Calculates and plots
}
