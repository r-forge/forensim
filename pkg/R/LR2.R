LR2<-function(Repliste,Tp,Td,Vp,Vd,xp,xd,theta,prDHet,prDHom,prC,freq){

#probbility of the evidence under H, function LRmix
	LRmix<-	function(Repliste,T,V,x,theta,prDHet,prDHom,prC,freq)
	{
		
		
		# function 1
		#    ProbEVid: function to calculate the probability Pr(U,T,V|H)
		ProbEvid <-	function(stain, freq, x=0, T = NULL, V = NULL, theta = 0 ) 
		{
			#require(genetics)
			if(!require(genetics)) stop("genetics package is required. Please install it.")
			if(!is.numeric(x) || is.na(x) || x < 0)
			{
				stop("'x' must be a positive integer")
			}
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
			
			if(length(theta)!=1)
			{
				stop("theta must be of length 1")
			}
			
			if (theta >= 1 || theta < 0) 
			{
				stop("'theta' must be a number in [0,1[" )
			}
			
			c <- length(stain)
			
			if(c != length(union(stain, stain)))
				stop("allele duplicates in 'stain'")
			
			# if (c != length(freq))
				# stop("'stain' and 'freq' must have the same length")
				
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
			# print(V)
			# if (sum(freq) == 1 && setequal(union(stain, allele.names(V)), stain) == FALSE)
				# stop("additional alleles in V (mixture contains all alleles of a locus)")
				
			# if (setequal(intersect(allele.names(T), stain), allele.names(T)) == FALSE)
				# stop("unknown alleles in 'T'")  
				
			if (theta >= 1 || theta < 0){stop("'theta' must be a number between 0 and 1, recommended 0.01 - 0.03")}
				
			# the known number of declared contributors to the mixture
			nT <- length(T)
			# the known number of people declared not to be contributors to the mixture 
			# (people who carry at least one allele from 'alleles')
			if (is.null(V) == FALSE) {nV <- length(V)}
			else{nV <- 0}

			
			   #theta correction: nV: number of non-contributors, 
			   #previsous version: only people with alleles shared with the known contributors are included, in this version, they are all included , this follows from curran et al 2007
			   # f <- 0
				# rep(1+,length(V))
				# for (i in 1:length(V)) 
				# {
					# print(sum(carrier(V[i], stain) == TRUE))#if the 
					# if (sum(carrier(V[i], stain) == TRUE) > 0)#if the ith non-contributor carries at least one allele that is in the sample, then add 1 to f
						# f <- f + 1#add 1 as one known non-contributor
				
				# print('f');print(f)
			
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
			vi <- rep(0, c)#vector of length c, the number of alleles in the sample: # of m alleles
			if (nV > 0) 
			{
				for (i in 1:c) 
				{
					vi[i] <- sum(allele.count(V, c(stain[i])))#known numbers of copies of Ai in V (sampling fromula <=> m=m+1)
				} 
			}
			# print(vi)
			# the known number of distinct stain carried by nT declared contributors
			t <- 0
			if (nT > 0)
				t <- length(allele.names(T))
			# the known number of distinct alleles which must be in U
			u<-c#u <- c - t
			# the known number of unknown contributors 
			nU <- x
			# the known number of alleles in U that can be of arbitrary type from C
			r <- 2 * nU #- u
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
			# if (nT > 0) 
			# {
				for (i in 1:c) {
					# if (sum(allele.count(T, stain[i])) == 0)  
						# ui[, i] <- ri[,i] + 1 ## because at least one copy from the T is there
					# else    
						ui[, i] <- ri[, i]
				}
			# }	   
			# else
				# ui <- ri + 1
			if (nU == 0)
				const <- factorial(2 * nU)
			else  
				const <- factorial(2 * nU) / prod((1 - theta) + ((2 * nT + 2 * nV):(2 * nT + 2 * nV + 2 * nU - 1)) * theta) 
			colnames(ui)<-stain
			results <- rep(0, nrow(ui))
			#print(length(results))
			loc1<-NULL
			for (d in 1:nrow(ui)) 
			{ 
				prod_p <- rep(0, c)
				for (i in 1:c) 
				{
					#if there are no copies of Ai (allele of type i) among the unknown contributors, then the probability pi is zero 
					#this is different from the original formula of Curran et al., here we define the conditional profile probability as the probability of the profile
					#under a certain hypothesis stating who gave the observed alleles, hence, Pr(stain="A"|U=0,V=0,T="A/A",H="suspect A/A gave the profile) would equal 1
					#rather then 2*p(A)*p(A) in the original formula
					if (ui[d, i] == 0){prod_p[i] <- 1}
					else{ 
					prod_p[i] <- prod((1 - theta) * freq[i] + ((ti[i] + vi[i]):(ti[i] + vi[i] + ui[d, i] -1)) * theta) 
					#if ui[i,j] is not null then print
					# loc1<-rbind(loc1,cbind(d,colnames(ui)[i],ui[d,i]))
					#loc1<-cbind(loc1,d,colnames(ui)[i])
					#print(class(ui[d,i]))
					}
				}
				
				results[d] <- prod(prod_p) / prod(factorial(ui[d, ]))
			}
			const*results

		}


		# auxiliary function 2
		# infoRep: function to calculate Pr(Ri|ui,T,V,H)
		

		infoRep <-function(R,Ggi,prDHet,prDHom,prC,freq)
		{
			#function to test wether the genotype is ho;ozygous
			homo<-function(gen)
			{
				if(gen[1]==gen[2]){TRUE}
				else{FALSE}
			}
			## the drop-in alleles set
			chi<-setdiff(R, intersect(R, as.numeric(Ggi)))
			# print('chi');print(chi)
			if(length(chi)==0)
			{	
				chifreq<-1-prC
			}
			#if chi is empty, then rho will take care of multiplying by (1-Prc), if I do it for chi
			else
			{
				chifreq<-rep(0,length(chi))#create a vector of length the number of contaminant alleles
				for(m in 1:length(chi))
				{
					#for every element of chi, write Pr(C)
					tmp1<-length(chi[m])#number of contaminant alleles
					#see curran et al pp49 for the justification of the factorial term
					chifreq[m]<-factorial(tmp1)*prod((prC^(tmp1))*freq[as.character(chi[m])])
				}
			}
			
			#dropout set
			i<-1
			#while loop to record the drop-outs per genotype (and so duplicates of the same dropped-out allele are allowed, though suhch scenario is unlikely in practice)
			rhovec<-vecdelta<-NULL
			while(i<=length(Ggi))
			{
				#which alleles have droppedout from genotype i
				delta2<-setdiff(as.numeric(Ggi[i:(i+1)]), intersect(R, as.numeric(Ggi[i:(i+1)])))
				# I use a temporary variable to record contamination on a genotype per genotype comparison basis
				chitemp<-setdiff(R, intersect(R, as.numeric(Ggi[i:(i+1)])))

				#which alleles haven't dropped-out
				# print('delta2');print(delta2)
				rho<-setdiff(R, union(as.numeric(chitemp), delta2))
				# print('rho');print(rho)
				##in case no drop-out is recorded, this is done here to avoid browsing Ggi twice
				if(length(rho)!=0)
				{
					
					##--- this corresponds to the rho vector is Curran's paper
					## we are comparing a gneotype (2 alleles) to the crime scene profile,if both aren't dropped out, record two, 
					##if it's a homozygote non-dropout
					if(homo(Ggi[i:(i+1)]))
					{ 
						deltarho<-prDHom
					}
					else#(!homo(Ggi[i:(i+1)]))
					{ 
						deltarho<-prDHet^length(rho)
					}
					##it's a heterozygote non-dropout
					rhovec<-c(rhovec,deltarho)
					delta<-1#
					# print('deltarho');print(deltarho)

					#at the end of the program: 1-rho
				}
				
				##if only one allele from the genotype has dropped-out, then it's either a  heterozygote or a homozygote, we have to treat both cases
				if(length(delta2)==1)
				{
					if(homo(Ggi[i:(i+1)])) {delta<-prDHom}
					if(!homo(Ggi[i:(i+1)])) {delta<-prDHet}
				}
				
				## in case it's a double het drop-out
				if(length(delta2)==2)
				{
					delta<-prDHet*prDHet
				}
				
				##record the drop-out
				vecdelta<-c(vecdelta,delta)
				## increment i
				i=i+2
			}
			
			 # print('rhovec');print(rhovec)
			 #in case all alleles have dropped-out
			if(is.null(rhovec)) rhofinal<-1
			else rhofinal<-1-prod(rhovec)#else, 1-prod[probability of dropout for the obsereved alleles]
			 # print(	prod(vecdelta)*(rhofinal)*chifreq)

			prod(vecdelta)*(rhofinal)*prod(chifreq)
			
		}

		
		
		Qdesign<-TRUE
		require('forensim')
		# the elements of U are detemined by the alleles across the replicates

		#the alleles/genotypes of the U set, there (n+r-1)!/(n-1)!r! possible ways of
		#choosing the genotypes for the unkown contributor(s), this is already implemented in forensim

		
		#-- if the Q designation is on
			if(Qdesign)
			{

				#get the set of unique alleles across the replicates
				
				alrep0<-unique(sort(unlist(Repliste)))
				# if(!is.null(alrep0))
				# {
					#Q: is the alleles that are observed in the population, but not in the samples
					Q<-setdiff(names(freq),alrep0)
					#update the list of possible alleles with Q
					alrep<-sort(unique(c(alrep0,Q)))
					# the number of alleles to sample from is given by the replicates length and the Q alleles
					C<-length(alrep0)
					# print(alrepQ)
					# print(1);print(Cmn(2*x,C ))
					d<-Cmn(2*x,C +length(Q))
					# print(2);print(d)
					#define set U
					Uset<-vector('list', length=d)
					#function comb gives the actual genotypic combination for the U set
					mat<-comb(2*x,C +length(Q))
					#now, we just have to rename the columns, accordingly to the allele names
					colnames(mat)<-alrep
					#now, we have to browse the matrix, and whenevr there is an allele, we replace
					# it with its name (taken form the corresponding column)
				# }
			}
			
			
			
		##---if the Q designation is shutdown
			else
			{
				#get the set of unique alleles across the replicates
				alrep<-unique(sort(unlist(Repliste)))
				if(length(alrep)!=0)
				{
					# print('ici')
					#Q: is the alleles that are observed in the population, but not in the samples
					# the number of alleles to sample from is given by the replicates length 
					C<-length(alrep)
					
					d<-Cmn(2*x,C)
					# print(d)
					#define set U
					Uset<-vector('list', length=d)
					#function comb gives the actual genotypic combination for the U set
					mat<-comb(2*x,C)
					#now, we just have to rename the columns, accordingly to the allele names
					colnames(mat)<-alrep
				}
			}
			#vector of frequencies for the unkown genotypes: Pr(Uj|T,V)
			# print('ici')
			# print(stain);print(freq);print(T);print(V);print(x);print(theta)
			# cat('alrep:',alrep,'\n')
			if(length(alrep)!=0)
			{
			a<-ProbEvid(stain=(alrep) ,freq=freq[as.character(alrep)],T=T,V=V,x=x,theta=theta)
			 # print(a);
			for(i in 1:d){Uset[[i]]<-rep(names(which(mat[i,]!=0)),mat[i,which(mat[i,]!=0)])}
			# print(Uset)

			#create the matrix of allele counts for the unkown contributors
			#whenevr there is an allele, we replace
			#it with its name (taken form the corresponding column)
			#as a reuslt, we have a list of the possible gneotypes for the unkown contributors
			# note that replicate information does not depend on people who have not contributed, 
			#so the V set 	does not interven in the calculation of the replicate probability
			# Gg: T union Uj
			Gg<-vector('list', length=d)
			#T could contain one or several contributors, obtain the list of alleles from T
			T2<-sapply(T,function(k) strsplit(k,'/'))
			# Gg2<-c(T,
			#first create the genotypes vector (refrence + unkown), that will be compared to 
			#the replicates
			if(length(T)!=0)#case where T is NULL and set to 0
			{
				for(i in 1:d)
				{
					Gg[[i]]<-unlist(c(T2,Uset[[i]]))#previsusly used union which did not account for duplicates
					# print('T2');print(T2);print('Uset');print(Uset[[i]])
				}
			}
			else
			{
				for(i in 1:d)
				{
					Gg[[i]]<-unlist(c(Uset[[i]]))#previsusly used union which did not account for duplicates
					# print('T2');print(T2);print('Uset');print(Uset[[i]])
				}
			
			}
			# print('Gg');print(Gg)
			resGeno<-rep(0,length=d)
			loc.delta<-NULL
			for(k in 1:length(Gg))#length(Gg)=d
			{
						
				# print('k');print(k)
				#for a given genotypic combination, calculate the product of the replicates probabilities
				Ggk<-as.numeric(Gg[[k]])#genotype k from the replicate
				#for a given Gg=(Ggk union Uk) combination, evaluate the probabilities for all the replicates 
				
				# print(Ggk)#;cat(freq)
				# print('prDHet');print(prDHet);print(prDHom);print('prDHom')
				# print('freq');	print(freq)
				resRep<-rep(0,length(Repliste))
				for(m in 1:length(Repliste))
				{
					
					# print('resrEp');print(infoRep8(Repliste[[m]],Ggk,prDHet,prDHom,prC,freq=freq))
					resRep[m]<-infoRep(Repliste[[m]],Ggk,prDHet,prDHom,prC,freq=freq)
					# print(Repliste[m])
					# print(resRep)
				}
				resRep<-prod(resRep)
				g<-(unique(Ggk))
				g2<-as.numeric(g)
				# cat('resRep\n',resRep,'\n')
				
				#if there are no unkwon contributors to the sample, then just return 1 x Pr(replicate)
				if(x==0)
				{	 
					resGeno[k]<-resRep
				}
				else
				{	#in case x!=0, then we must sum the Pr(Ri|Ui) for each possible Ui
					
					resGeno[k]<-a[k]*resRep
				}
			}
			 sum(resGeno)
			  # write.table(Uset,file='Uset.txt')
			}
			else {
			# print('oui')
			return(0)}

		 
	}




num<-LRmix(Repliste,T=Tp,V=Vp,x=xp,theta,prDHet,prDHom,prC,freq)
deno<-LRmix(Repliste,T=Td,V=Vd,x=xd,theta,prDHet,prDHom,prC,freq)
list('num'=num,'deno'=deno,'LR'=num/deno)
}