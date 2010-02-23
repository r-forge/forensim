#Hinda Haned, November 2009

tabSPH<-function(x,y=NULL,geno,tabcsv,byloc=FALSE,s=40)
{
	tab.loc<-recordHeights(x,y,geno,tabcsv,byloc)
	
	
	if(byloc)
	{
		m<-names(geno)
		perloc<-vector('list',length(m))
		names(perloc)<-m
		for(M in names(tab.loc))
		{
			
			#the case where no geneotypes are selected, in case of allele sharing yields a null data.frame, is ignored
			#only non empty data frames are of interest here
			#print(tab.loc[[M]])
			if(length(tab.loc[[M]])!=0)
			{
				#if both peak heights are null, we are not really inetrested
				if(!setequal(c(tab.loc[[M]][,"Height1"]!=0,tab.loc[[M]][,"Height2"]!=0),c(0,0)))
				{
					tmp1<-NULL
					tmp2<-NULL
					for(i in tab.loc[[M]][,"Height2"])
					{
						if(i>s) tmp1<-c(tmp1,1)# no drop-out, peak height >40
						else{ tmp1<-c(tmp1,0)}#allele dropped out
					}	
					tab1<-cbind.data.frame("D"=tmp1,"SPH"=tab.loc[[M]][,"Height1"])
					#if the height of the first allele is not null then we score 1, the height of allele 2 is stored
					for(i in tab.loc[[M]][,"Height1"])
					{
						if(i>s) tmp2<-c(tmp2,1)# no drop-out, peak height >40
						else{ tmp2<-c(tmp2,0)}#allele dropped out
					}	
					tab2<-cbind.data.frame("D"=tmp2,"SPH"=tab.loc[[M]][,"Height2"])
					res<-rbind(tab1,tab2)
					#we want to model the dropout probability, so the first column should be equal to one in case of dropout
					res[,1]<-1-res[,1]
					perloc[[M]]<-res
				}
			}
			#else{perloc[[M]]<-cbind.data.frame('d'=0,'h'=0)}
		}
		perloc

	}
	#if data overall loci, things get simpler: tab.loc is a matrix with two columns Height1 & Height2
	else
	{
		
		#if the height of the second allele is not null then we score 1, the height of allele 1 is stored
		tmp1<-NULL
		tmp2<-NULL
		for(i in tab.loc[,"Height2"])
		{
			if(i>s) tmp1<-c(tmp1,1)# no drop-out, peak height >40
			else{ tmp1<-c(tmp1,0)}#allele dropped out
		}	
		tab1<-cbind.data.frame("D"=tmp1,"SPH"=tab.loc[,"Height1"])
		#if the height of the first allele is not null then we score 1, the height of allele 2 is stored
		for(i in tab.loc[,"Height1"])
		{
			if(i>s) tmp2<-c(tmp2,1)# no drop-out, peak height >40
			else{ tmp2<-c(tmp2,0)}#allele dropped out
		}	
		tab2<-cbind.data.frame("D"=tmp2,"SPH"=tab.loc[,"Height2"])
		res<-rbind(tab1,tab2)
		#we want to model the dropout probability, so the first column should be equal to one in case of dropout
		res[,1]<-1-res[,1]
		res
		
	}
}
		
		
