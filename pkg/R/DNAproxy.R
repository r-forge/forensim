#Hinda Haned, November 2009
DNAproxy<-function(tab,x)
{
	#if tab is generated form DNA mixtures, then the number of columns must be 6
	#this a temporary way to test for mixtures, to be improved !
	if(length(tab[[1]])==6)
	{
		mix<-TRUE
	}
	else{mix<-FALSE}
	DNAproxy.loc<-function(tabloc,x)
	{
		I<-intersect(tabloc$Expected,tabloc$Observed)
		if(length(I)!=0)
		{
			n<-rep(0,length(I))
			names(n)<-I
			H<-tabloc[,"Heights"]
			names(H)<-tabloc[,"Expected"]
			h<-n; s<-n
			if(mix)
			{
				#x1  are the non dropped out alleles, belonging to the first individual
				x1<-tabloc[tabloc$Expected %in% I,x]
				#x2  are the non dropped out alleles, belonging to the second individual
				x2<-tabloc[tabloc$Expected %in% I,names(tabloc)[2:3][which(!(names(tabloc)[2:3]) %in% x)]]
				for(i in 1:length(I))
				{	
					h[I[i]]<-H[I[i]]#h stores the peak heights of the non dropped out alleles
					
					if(x1[i]==0)
					{	n[I[i]]<-0#n records the number of hetero and homo alleles
						s[I[i]]<-0#single contributed alleles, s give the indexes of the elements to be taken into account in the calculation
					}
					if(x1[i]!=0)
					{ 
						if(x1[i]==1 & x2[i]==0){n[I[i]]<-1; s[I[i]]<-1}
						if(x1[i]==1 & x2[i]>1){n[I[i]]<-1; s[I[i]]<-0}
						if(x1[i]==2 & x2[i]==0){n[I[i]]<-2; s[I[i]]<-1}
						if(x1[i]==1 & x2[i]==1){n[I[i]]<-1; s[I[i]]<-0}
						if(x1[i]==2 & x2[i]>1){n[I[i]]<-2; s[I[i]]<-0}
					}
				}
				cbind.data.frame('deno'=sum(n*s),'num'=sum(h*s))
			}
			#the procedure is straightforward if the data is composed of single contributor stains
			else
			{
				#all alleles are taken into account, we don't need the s vector
				#n records the coefficients : 1 for alleles carried by a heterozygote and 2 for alleles carried by a homozygote
				#h records the peak heights of the observed alleles
				for(i in 1:length(I))
				{
					#if two alleles are expected, then all elements of n are equal to 1
					if(length(tabloc$Expected)==2) {n[I[i]]<-1}
					#if a single allele is expected then  length(I) is one at most, and n[I[i]]<-2
					if(length(tabloc$Expected)==1 ){n[I[i]]<- 2}
					h[I[i]]<-H[I[i]]
				}
				#sum(n): n=2 if allele carried by a homozygote, n=1 if allele carried by a heterozygote; sum(h): sum of the peak heights observed
				cbind.data.frame('deno'=sum(n),'num'=sum(h))
			}
			
		}

		#if all alleles have dropped out, the locus can not be taken into account for the DNA proxies calculations
		else#if all alleles have dropped out, the locus can not be taken into account for the DNA proxies calculations
		{
			cbind.data.frame('deno'=0,'num'=0)
		}
	}

	#tabo<-org(x,y,geno,tab)
	#DNAproxy.loc(tab[[1]],'c1')
	ind<-sum(unlist(sapply(tab,DNAproxy.loc,x)[2,]))/sum(unlist(sapply(tab,DNAproxy.loc,x)[1,]))
	if(is.na(ind)){ind<-0}
	ind
}



	



