\name{dropDB}
\alias{dropDB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates LR allowing for drop-out and drop-in
}
\description{
This function wraps up David Balding's code.
}
\usage{
dropDB(dataBase, DNA, dropHetero = 0.05, alpha = 0.5, 
dropIn = 0.05, rel = c(0, 0)/4, maxUnknown = 1,
 adj = 0, fst = 0)
}

\arguments{
  \item{dataBase}{
An R data frame with header 'marker allele count'
with columns giving the names of the marker and the allele and the count
(the number of occurences of the allele) in the database.
}
  \item{DNA}{
An R data frame with header 'MarkerSource A1 A2 A3 A4'
There are four lines for each marker corresponding to alleles
for the (i) the mixture (ii) known contributor (if there is one) (iii) the defendant and
(iv) drop-out alleles. The names and order of the markers should be as for dataBase-} 

  \item{dropHetero}{
Probability heterozygous drop-out
}
  \item{alpha}{
Balding and Buckleton's alpha-parameter defining relation between
homozygous and heterozygous drop out probabilities.
}
  \item{dropIn}{
Drop-in probability.
}
  \item{rel}{
A vector of length 2 giving the IBD=1 and IBD=2 probabilities
and so c(0,0) corrresponds to unrelated contributors.
}
  \item{maxUnknown}{
The maximum number of unknown contributors under the defence hypothesis.
Possible values 1 or 2.
}
  \item{adj}{
Parameter to adjust for sampling adjustment. Balding uses 2.
}
  \item{fst}{
Parameter to adjust for coancestry. Balding recommends values between 0.01 and 0.05.
}
}

\value{
Marker name, LR and numerator and denominator of LR.
}
\references{
Balding DJ, Buckleton J, Interpreting low template DNA profiles,
Forensic Science International: Genetics, 4: 1-10, 2009, doi: 10.1016/j.fsigen.2009.03.003.
http://www.zebfontaine.eclipse.co.uk/djb.htm
}
\author{
Thore Egeland based on David Balding's code.
}





\examples{
#First output column of Table 1 in Balding and Buckleton:
data(Bates.Database);data(Bates.DNA)
dropDB(Bates.Database,Bates.DNA,dropHetero=0.05,alpha=0.5,dropIn=0.00,rel=c(0,0)/4,
maxUnknown=1,adj=2,fst=0.02)
#First output column of Table 3 in Balding and Buckleton:
data(CaseY.Database);data(CaseY.DNA)
dropDB(CaseY.Database,CaseY.DNA,dropHetero=0.5,alpha=0.5,dropIn=0.05,rel=c(0,0)/4,
maxUnknown=1,adj=2,fst=0.02)
}
