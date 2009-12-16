#  Hinda Haned, June 9 2009

Pevid2 <-
function(stain, freq, x=0, T = NULL, V = NULL, theta = 0 ) 
{
    #require(genetics)
	if(!require(genetics)) stop("genetics package is required. Please install it.")
	if(!is.numeric(x) || is.na(x) || x < 0)
	{
        stop("'x' must be a positive integer")
    }
 	x0 <- round(x)
    if (x != x0) 
		warning("'x' has been rounded to integer")
    if (!is.vector(stain)) 
        stop("'stain' must be a vector")
		
	if (round(sum(freq)) > 1) 
	{
		stop("sum of allele frequencies must not exceed 1")
	}
		
	if (!is.numeric(freq) || any(is.na(freq)) || any(freq < 0) || any(freq > 1))
	{
		stop("stain frequencies in  'freq' must be numbers strictly between 0 and 1")
	}
	
	if (!is.vector(freq)) 
	{
        stop("'freq' must be a vector")
	}
	
   	if (theta >= 1 || theta < 0) 
	{
		stop("'theta' must be a number in [0,1[" )
	}
    
	c <- length(stain)
	
    if(c != length(union(stain, stain)))
        stop("allele duplicates in 'stain'")
	
    if (c != length(freq))
        stop("'stain' and 'freq' must have the same length")
		
    if (!is.numeric(freq) || any(is.na(freq)) || any(freq <= 0) || any(freq >= 1))
        stop("all entries of 'freq' must be numbers between 0 and 1")
		
  	if (round(sum(freq)) > 1) 
	{
			stop("sum of allele frequencies must not exceed 1")
	}
		
    if (is.null(T) == FALSE)
	{
        if (is.genotype(T) == FALSE)
            T <- as.genotype(T)
        if (any(is.na(T)))
            stop("entries of T are not correct")
    }      
    if (is.null(V) == FALSE)
	{
        if (is.genotype(V) == FALSE)
            V <- as.genotype(V)
        if (any(is.na(V)))
            stop("entries of V are not correct")
    }
    if (sum(freq) == 1 && setequal(union(stain, allele.names(V)), stain) == FALSE)
        stop("additional alleles in V (mixture contains all alleles of a locus)")
		
    if (setequal(intersect(allele.names(T), stain), allele.names(T)) == FALSE)
        stop("unknown alleles in 'T'")  
		
    if (theta >= 1 || theta < 0)
	{
        stop("'theta' must be a number between 0 and 1, recommended 0.01 - 0.03") 
    }
	
    # the known number of declared contributors to the mixture
    nT <- length(T)
    # the known number of people declared not to be contributors to the mixture 
    # (people who carry at least one allele from 'alleles')
    if (is.null(V) == FALSE) 
	{
        f <- 0
        for (i in 1:length(V)) 
		{
			if (sum(carrier(V[i], stain) == TRUE) > 0)
				f <- f + 1
		}
		nV <- f
    }
    else
        nV <- 0
    # the known number of heterozygous declared contributors (in T)
    if (nT > 0)
        hT <- sum(heterozygote(T) == TRUE)  
    else 
        hT <- 0  
    # the known number of heterozygous declared non-contributors (in V)
    if (nV > 0)
        hV <- sum(heterozygote(V) == TRUE)  
    else 
        hV <- 0
    # the known number of copies of allele 'stain[i]' in T
    ti <- rep(0, c)
    if (nT > 0)
	{
        for (i in 1:c)
		{
            ti[i] <- sum(allele.count(T, stain[i]))
        }
    }
    # the known number of copies of allele 'stain[i]' in V
    vi <- rep(0, c)
    if (nV > 0) 
	{
        for (i in 1:c) 
		{
            vi[i] <- sum(allele.count(V, stain[i]))
        } 
    }
    # the known number of distinct stain carried by nT declared contributors
    t <- 0
    if (nT > 0)
        t <- length(allele.names(T))
    # the known number of distinct alleles which must be in U
    u <- c - t
    # the known number of unknown contributors 
    nU <- x
    # the known number of contributors
    nC <- nT + nU
    min.nU <- trunc(u / 2) + u %% 2
    min.nC <- trunc(c / 2) + c %% 2 
    if (nU < min.nU || nC < min.nC)
        stop("not enough contributors")
    # the set of distinct alleles which must be in U:  U0 = C \ T
    if (nT > 0)
        U0 <- setdiff(stain, allele.names(T))
    else
        U0 <- stain
    # the known number of alleles in U that can be of arbitrary type from C
    r <- 2 * nU - u
    # ri[,i] ... the unknown number of copies of allele 'alleles[i]' among 
    # the r unconstrained alleles in U
    if (r>0)
	{ 
        # all possibilities of the vector (r_1, r_2, ... , r_c), where  
        # sum(ri)=r and ri>=0, ri<=r
        ri <- comb(r,c)# previously, the  perm function, comb calls a C function, recurs, see the src directory in package source
    }
    else
	{
        ri <- matrix(rep(0, c), nrow = 1, byrow = TRUE)
    }
	# ui[,i] ... the unknown number of copies of allele 'alleles[i]' in U, 
    # sum(ui) = 2*nU 
    ui <- matrix(0, nrow = nrow(ri), ncol = c)
    if (nT > 0) 
	{
        for (i in 1:c) {
            if (sum(allele.count(T, stain[i])) == 0)  
                ui[, i] <- ri[,i] + 1
            else    
                ui[, i] <- ri[, i]
        }
    }	   
    else
        ui <- ri + 1
    if (nU == 0)
        const <- factorial(2 * nU)
    else  
        const <- factorial(2 * nU) / prod((1 - theta) + ((2 * nT + 2 * nV):(2 * nT + 2 * nV + 2 * nU - 1)) * theta) 
    results <- rep(0, nrow(ui))
    for (d in 1:nrow(ui)) 
	{ 
        prod_p <- rep(0, c)
        for (i in 1:c) 
		{
            #if there are no copies of Ai (allele of type i) among the unknown contributors, then the probability pi is zero 
			#this is different from the original formula of Curran et al., here we define the conditional profile probability as the probability of the profile
			#under a certain hypothesis stating who gave the observed alleles, hence, Pr(stain="A"|U=0,V=0,T="A/A",H="suspect A/A gave the profile) would equal one
			#rather then 2*p(A)*p(A) in the original formula
 			if (ui[d, i] == 0){prod_p[i] <- 1}
			else{ prod_p[i] <- prod((1 - theta) * freq[i] + ((ti[i] + vi[i]):(ti[i] + vi[i] + ui[d, i] -1)) * theta) }
        }
		
        results[d] <- prod(prod_p) / prod(factorial(ui[d, ]))
    }
    return(const * sum(results))
}

