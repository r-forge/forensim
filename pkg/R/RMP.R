#  Hinda Haned, June 9 2009
# Random Match Probability (RMP)

RMP <-
function(suspect=NULL, filename=NULL, freq, k = c(1, 0, 0), theta = 0,refpop=NULL)
{
	if(!require(genetics)) stop("genetics package is required. Please install it.")
	if(!is.tabfreq(freq) & !is.vector( freq))
	{
        stop("'freq' is either a tabfreq objcet (multlocus case) or a vector (single locus case)")
    }
   
		
    if (!is.numeric(k) || any(is.na(k)) || any(k < 0) || any(k > 1)) 
        stop("all entries of 'k' must be numbers from the interval [0, 1]")
		
    if (length(k) != 3)
        stop("'k' should be a vector containing three kinship coefficients")

	if (sum(k) != 1) 
        stop("kinship coefficients mut some to 1")

	if(!is.null(filename) & !is.null(suspect))
		warning('only data in suspect argument will be used')
		tab <- suspect
		
	if(is.null(filename) & !is.null(suspect) )
		tab <-suspect
	
	if(is.null(suspect) & !is.null(filename))
	
		tab <- read.table(filename,header=TRUE)
	
	if(!is.matrix(tab) & !is.data.frame(tab))
		stop(cat("data in filename", " must be a matrix or a data frame"))
	nc <- ncol(tab)
	
	if(nc!=2)
		stop(cat("data in must be a matrix or a data frame of dimension 2 x number of loci " ))
	
	locus <- as.character(tab[,1])
	geno <- as.character(tab[,2])
	#print(tab)
	names(geno) <- locus
	#print(geno)
	#geno2 <-as.genotype(geno)
	freq.loc <- freq$which.loc
	if(!any(locus %in% freq.loc))
	{
		stop("There are uknown loci in the DNA stain, please check your data")
	}
	pop <-freq$pop.names
	
	if(!is.null(pop) & length(pop)==1)
	{
		AF <- freq$tab[[pop]]
	}
	
	if(length(pop)>1)
	{
		AF <- freq$tab[[refpop]]
	}
	
	if(is.null(refpop))
		AF <- freq$tab[locus]
	
	#getting the corresponding allele frequencies pi
	pi <- NULL
	
	for(l in locus)
	{
		pi <- rbind(pi,(AF[[l]][strsplit(geno[l],'/')[[1]]]))
	}	
	
	colnames(pi) <- NULL
	
	
	#heterozygote state for each gneotype at each locus
	hetstate<-NULL
	for(i in geno)
	{	
		loc <- strsplit(i,'/')[[1]]
		if(loc[1]==loc[2]){hetstate <-c(hetstate,FALSE)}
		else{hetstate<-c(hetstate,TRUE)}
	}
	
	
	res <- vector('list',length(locus))
	
    for (h in 1:length(locus))#or in 1:hetstate, as hetstate is of length nlcous
	{
        #print('hstate');print(hetstate[h])
		if(!hetstate[h])#homo case
		{
			res[[h]] <- k[1]*(2*theta + (1 - theta)*pi[h,1])* 
			(3*theta + (1 - theta)*pi[h,1])/((1+theta)*(1+2*theta)) +
			k[2]*(2*theta + (1 - theta)*pi[h,1])/(1 + theta) +
			k[3]
		}
		
        else
		{
			res[[h]] <- k[1]*2*(theta + (1 - theta)*pi[h, 1])*
			(theta + (1 - theta)*pi[h, 2])/((1+theta)*(1+2*theta))+ 
			k[2]*(2*theta + (1 - theta)*(pi[h, 1] + pi[h,2 ]))/(2*(1 + theta)) + 
			k[3]	
		}
	}
	names(res) <-locus
    final <-vector('list',2)
	names(final) <- c('RMP.loc','RMP')
	final$RMP.loc <-signif(unlist(res),2)
	final$RMP <- signif(prod(unlist(res)),2)
	return(final)
	
}

