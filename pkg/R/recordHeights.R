#Hinda Haned, November 2009
#simply displays in a matrix of two columns the peak heights observed for each individual, when the considered genotype is heterozygous 
#two heterozyhotes with non shared alleles are also taken into account 
recordHeights<-function(x,y=NULL,geno,tabcsv,byloc=FALSE)
{
    tabo<-recordDrop(x,y,geno,tabcsv)
	#displays the dropouts events, the genotypes of both contributors and the corresponding peak heights
	#DNAproxy estimations
	#if tabo is generated form DNA mixtures, then the number of columns must be 6
	#this a temporary way to test for mixtures, to be improved !
	if(length(tabo[[1]])==6)
	{
		mix<-TRUE
	}
	else{mix<-FALSE}
	resloc<-vector('list',length(tabo))#list of matrices where peak height1 and peak height2, in rows for each individual
	names(resloc)<-names(tabo)
	#in case of mixed DNA stains
	if(mix)
	{
		for(i in 1:length(tabo))
	    {
	        a1<-NULL;a2<-NULL
			mat<-tabo[[i]]
			#only non shared alleles carried by heterozygote individuals are included in the calculation
			cont1<-t(matrix(mat[mat[,paste("c",x,sep="")]!=0 & mat[,paste("c",y,sep="")]==0,"Heights"]))
			cont2<-t(matrix(mat[mat[,paste("c",y,sep="")]!=0 & mat[,paste("c",x,sep="")]==0,"Heights"]))
			if(length(cont1)!=0 & ncol(cont1)==2)#in the non homozygote case
	        {
	            a1<-cbind.data.frame(cont1[,1],cont1[,2])
				colnames(a1)<-c("Height1","Height2")
	        }
	        if(length(cont2)!=0 & ncol(cont2)==2)
	        {
	            a2<-cbind.data.frame(cont2[,1],cont2[,2])
				colnames(a2)<-c("Height1","Height2")
	        }
	        resloc[[i]]<-rbind.data.frame(a1,a2)
			
		}
                
    }
	#in case of single-contributor stains
    else
	{
		for(i in 1:length(tabo))
	    {
	        
			a1<-NULL;a2<-NULL
			mat<-tabo[[i]][,'Expected']
			#if length(mat)==2, then we are delain with a hetrozygote
			if(length(mat)==2)
			{
				if(tabo[[i]][1,'Heights']!=0 & tabo[[i]][2,'Heights']!=0)#the case where all peak heights are null is skipped
				{
					a1<-c(a1,tabo[[i]][1,'Heights'])
					a2<-c(a2,tabo[[i]][2,'Heights'])
					resloc[[i]]<-cbind.data.frame('Height1'=a1,'Height2'=a2)
				}
			}
		}
	    
	}
	#if results are to be displayed per locus, no changes
	if(byloc){resloc}
	#otherwise, concatenate the results to obtain a single matrix with two columns: Height1 & Hieght2
	else
	{
		a1<-NULL
		a2<-NULL
	    for(i in 1:length(resloc))
	    {
	        if(length(resloc[[i]])!=0)
			{	
				a1<-rbind.data.frame(a1,matrix(resloc[[i]][,'Height1']))
				a2<-rbind.data.frame(a2,matrix(resloc[[i]][,'Height2']))
				               
			}
	   
	    	
		}
		matres<-cbind.data.frame(a1,a2)
		colnames(matres)<-c('Height1','Height2')
		#removing null rows (null peak heights for both alleles)
		selec<-which(matres[,1]==0 & matres[,2]==0)
		if(length(selec)!=0){ matres<-matres[-which(matres[,1]==0 & matres[,2]==0),]}
		matres
	}
}
