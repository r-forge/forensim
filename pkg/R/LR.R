# Hinda Haned, June 2 2009
#evaluation of forensic DNA stains through likelihood ratios
LR <-
function(stain, freq, xp=0,xd=0, Tp = NULL, Vp = NULL, Td=NULL,Vd=NULL, theta = 0 )
{
	LR1 <- Pevid2(stain=stain, freq=freq, x=xp, T=Tp, V=Vp, theta = 0 )
	LR2 <- Pevid2(stain=stain, freq=freq, x=xd, T=Td, V=Vd,theta = 0 )
	LR <- LR1/LR2
	return(LR)
  
}

#__________________________________________________________#
#readfile version, coming soon
# lik.ratio2 <-
# function(filename)
# {
	# fopen <- read.csv('filename',h=TRUE,sep="\t")
	# "xls2csv"
	# tmp <- colnames(fopen)
	# nloc <- nrow(fopen)#number of loci
	
	# LR1 <- Pevid2(stain=stain, freq=freq, x=xp, T = T1, V = V1, theta = 0 )
	# LR2 <- Pevid2(stain=stain, freq=freq, x=xd, T=T2, V=V2,theta = 0 )
	# LR <- LR1/LR2
	# return(LR)
  
# }