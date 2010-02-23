#Hinda Haned, November 2009
tabDNAproxy<-function(x,y=NULL,geno,tabcsv)
{
	tabo<-recordDrop(x,y,geno,tabcsv)
	#DNA proxy estimations
	#if tabo is generated form DNA mixtures, then the number of columns must be 6
	#this a temporary way to test for mixtures, to be improved !
	if(length(tabo[[1]])==6)
	{
		mix<-TRUE
	}
	else{mix<-FALSE}
	
	m<-names(geno)
	perloc<-vector('list',length(m))
	names(perloc)<-m
	
	#if the analyzed stains are 2-person mixtures
	if(mix)
	{
		ind1<-DNAproxy(tabo,paste("c",x,sep=""))
		ind2<-DNAproxy(tabo,paste("c",y,sep=""))
		
		for(M in names(tabo))#for each locus
		{
			Dloc<-NULL
			Hestim<-NULL
			for(h in 1:nrow(tabo[[M]]))
			{
				Dloc<-c(Dloc,tabo[[M]][h,"D"])
				Hestim<-c(Hestim,sum(tabo[[M]][h,2:3]*c(ind1,ind2)))
			}
			perloc[[M]]<-cbind.data.frame(Dloc,Hestim)
		}
	}
	#single contributors stains
	else
	{
		ind<-DNAproxy(tabo)
		for(M in names(tabo))#for each locus
		{
			Dloc<-NULL
			Hestim<-NULL
			for(h in 1:nrow(tabo[[M]]))
			{
				#print(tabo[[M]])
				Dloc<-c(Dloc,tabo[[M]][h,"D"])
				#hetrozygote state	
				if(length(tabo[[M]]$Expected)==1){ state <-2} # equals two if the genotype is homozygous
				if(length(tabo[[M]]$Expected)==2){ state<-1} #one otherwise
				Hestim<-c(Hestim,sum(state*ind))
			}
			perloc[[M]]<-cbind.data.frame(Dloc,Hestim)
		}
	}
	perloc
}
		
		
