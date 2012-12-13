#likelihood ratios new version allowing empty samples
#probbility of the evidence under H, function LRmix

# first define the  function which allows computing the likelihood of a given hypothesis

likEvid<-function(Repliste,T,V,x,theta,prDHet,prDHom,prC,freq)
{	
	
	infoRep <- function(R,Ggi,prDHet,prDHom,prC,freq)
	{
		# print('R');print(R);
		R<-as.numeric(R)
		# print(R)
		# print(prDHet);print(prDHom);print(prC);print(freq)
		# print('Ggi');print(as.numeric(Ggi)); print(Ggi)
		Ggi<-as.numeric(Ggi)
		# print(Ggi)
		res <-0
		allele<-names(freq)
		freq0<-as.numeric(freq)
		# print(R);print(Ggi);print(allele);print(freq0)

			infoRepC<-function(R,Ggi,prDHet,prDHom,prC,allele,freq0,res)
			{
				lenR<-length(R)
				lenGgi<-length(Ggi)
				lenFreq<-length(freq)
				lenHom<-length(prDHom)
				lenHet<-length(prDHet)
				
				.C('infoRepC',as.double(R),as.integer(lenR),as.double(Ggi),as.integer(lenGgi),as.double(prDHet),as.integer(lenHet),as.double(prDHom),as.integer(lenHom),as.double(prC),as.character(allele),as.double(freq0),as.integer(lenFreq),res=as.double(res),PACKAGE='forensim')
			}
			loc<-infoRepC(R,Ggi,prDHet,prDHom,prC,allele,freq0,res)
			# print(loc$res)
			loc$res
	}


	Pevid6 <- function(stain, freq, x, T, V, theta) 
	{
		# require('forensim')
		nbligne<-Cmn(2*x,length(stain))
		# factorial(n + m - 1)/(factorial(n - 1) * factorial(m))

		#pre-allocation of a vector of length stain, that will contain the output from the C function
		foo<-rep(0,nbligne)
		#essaye de la faire independemment?
		PevidC<-function(stain, freq, x, T, V,theta,foo)
		{
			lenStain<-length(stain)
			lenFreq<-length(freq)
			lenT<-length(T)/2
			lenV<-length(V)/2

			.C('PevidC',as.double(stain),as.integer(lenStain),as.double(freq), as.integer(lenFreq),as.integer(x),as.double(T),as.integer(lenT),as.double(V),as.integer(lenV),as.double(theta),as.integer(nbligne),foo=as.double(foo), PACKAGE='forensim')
			
		}

		loc<-PevidC(stain, freq, x, T, V,theta,foo)
		# print(loc$res)
	 loc$foo
	}


	
	require('forensim')
	# the elements of U are detemined by the alleles across the replicates
	#-- if the Q designation is on
		
	#get the set of unique alleles across the replicates
	alrep<-unlist(Repliste)
	if(length(alrep)==0)
	{
		Q<-names(freq)
	}
	else{
	#if we see at least one allele among all replicates
		alrep2<-unique(sort(alrep))#sort on empty set may generate warning

		#Q: is the alleles that are observed in the population, but not in the sample
		Q<-unique(c(names(freq),alrep2))
		#now, we have to browse the matrix, and whenevr there is an allele, we replace
		# it with its name (taken form the corresponding column)
	}

	# print(Q)
	d<-Cmn(2*x,length(Q))
	#define set U
	Uset<-vector('list', length=d)
	#function comb gives the actual genotypic combination for the U set
	mat<-comb(2*x, length(Q))
	#now, we just have to rename the columns, accordingly to the allele names
	colnames(mat)<-Q
	#now, we have to browse the matrix, and whenevr there is an allele, we replace
	# it with its name (taken form the corresponding column)
# }
	
##---if the Q designation is shutdown
	
	#vector of frequencies for the unkown genotypes: Pr(Uj|T,V)
	rare<-Q[which(!Q %in% names(freq),arr.ind=TRUE)]
	if(length(rare)!=0){
	freq[as.character(rare)]<-1/(2*2085)
	print(paste('WARNING: allele',rare,'has been added to the database with frequency',1/(2*2085),sep=' '))}
	# print(rare)
	a<-Pevid6(stain=Q ,freq=freq[as.character(Q)],T=T,V=V,x=x,theta=theta)
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
	# T2<-sapply(T,function(k) strsplit(k,'/'))
	# Gg2<-c(T,
	#first create the genotypes vector (refrence + unkown), that will be compared to 
	#the replicates
	if(length(T)!=1)#case where T is NULL and set to 0
	{
		for(i in 1:d)
		{
			Gg[[i]]<-unlist(c(T,Uset[[i]]))#previsusly used union which did not account for duplicates
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
		Ggk<-Gg[[k]]#genotype k from the replicate
		u<-Uset[[k]]#kth possible unkown genotype: Uk
		#for a given Gg=(Ggk union Uk) combination, evaluate the probabilities for all the replicates 
		# cat(Ggk,'\n')#;cat(freq)
		 # print('Ggk');print(Ggk)
		# print('prDHet');print(prDHet);print(prDHom);print('prDHom')
		# print('freq');	print(freq)
		resRep<-rep(0,length(Repliste))
		for(m in 1:length(Repliste))
		{
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


LR<-function(Repliste,Tp,Td,Vp,Vd,xp,xd,theta,prDHet,prDHom,prC,freq){

num<-likEvid(Repliste,T=Tp,V=Vp,x=xp,theta,prDHet,prDHom,prC,freq)
deno<-likEvid(Repliste,T=Td,V=Vd,x=xd,theta,prDHet,prDHom,prC,freq)


list('num'=num,'deno'=deno,'LR'=num/deno)}
