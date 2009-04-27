# Hinda Haned, April 2009
# haned@biomserv.univ-lyon1.fr

PE <-
function(mix, freq, refpop=NULL, theta=0,byloc=FALSE,digits=2)
{
	options(digits=digits)
	popinfo <-mix@popinfo
	loc <- mix@which.loc
	pe.loc<-matrix(0,nrow=1,ncol=length(loc))
	colnames(pe.loc) <- loc
	rownames(pe.loc) <- "PE_l"
	if(is.null(popinfo))
	{
		af <- findfreq(mix,freq)
		print(af)
	}
	else
	{
		af <- findfreq(mix,freq, refpop)[[refpop]]
		
	}
	
	if(theta==0)
	{
		for( l in  loc)
		{
			pe.loc[1,l]<-1-sum((af[[l]]*af[[l]]))
		}
		
		
		if(byloc)
		{
			return(pe.loc)
		}
		
		else
		{
			PE <- 1-prod(1-pe.loc)
			names(PE) <- "PE"
			return(PE)
		}
	}
	if(theta!=0)
	{
		if (theta >= 1 || theta < 0) 
		{
			stop("'theta' must be a number in [0,1[" )
		}
		
		for( l in  loc)
		{
			pe.loc[1,l]<-1-sum((af[[l]]*af[[l]]))-theta*sum(af[[l]])*(1-sum(af[[l]]))
		}
		
		
		if(byloc)
		{
			return(pe.loc)
		}
		
		else
		{
			PE <- 1-prod(1-pe.loc)
			names(PE) <- "PE"
			return(PE)
		}
	}
}








