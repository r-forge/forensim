#Hinda Haned, November 2009
recordDrop<-function(x,y=NULL,geno, tabcsv,s=40)
{
	#the genetics package is required in order to perform the allele counts
	if(!require(genetics)) stop("genetics package is required. Please install it.")
	loci<-names(geno)
	loci2<-tabcsv$Marker
	if(!all(loci %in% loci2))
	{	
		stop("different loci names in geno and tabcsv, loci names are case sensitive")
	}
	else
	{
		#initialize the variables
		res<-vector('list',length(loci))
		names(res)<-loci
	}
	#if data is a mixture, then y is  null
	if(is.null(y))
	{
		for(k in 1:length(loci))
		{
			#for each locus, store the allele counts for each person
			mat<-as.matrix(allele.count(genotype(as.character(as.matrix(geno[c(x,y),k])))))
			if(loci[k]!="AMEL" & loci[k]!="Amel")
			{
				if(ncol(mat)!=1)
				{
					tmp<-sort(as.numeric(colnames(mat)))
					tmp<-as.character(tmp)
					mat3<-t(as.matrix(mat[,tmp]))
				}
				if(ncol(mat)==1)
				{

					tmp<-sort(as.numeric(colnames(t(mat))))
					tmp<-as.character(tmp)
					mat3<-as.matrix(t(mat)[,tmp])

				}
			}
			if(loci[k]=="AMEL" || loci[k]=="Amel")
			{
				if(ncol(mat)==1)
				{
					
					tmp<-sort((colnames(t(mat))))
					tmp<-as.character(tmp)
					mat3<-as.matrix(t(mat)[,tmp])

				}
				if(ncol(mat)!=1)
				{
					tmp<-sort(colnames(mat))
					mat3<-t(as.matrix(mat[,tmp]))
				}
			}

			#if the data is not a mixture, only two alleles and two peak heights are expected (at most)
			#each element of res will store a  matrix of the results for allele counts expected and observed alleles and the peak heights
			res[[k]]<-cbind.data.frame(t(t(rownames(mat3))),t(tabcsv[tabcsv$Marker==loci[k],c("Allele1","Allele2")][1:nrow(mat)]),t(tabcsv[tabcsv$Marker==loci[k],c("Height1","Height2")][1:nrow(mat)]))
			#res[[k]]<-cbind.data.frame(rownames(restmp),restmp)
			colnames(res[[k]])<-c("Expected","Observed","Heights")
			rownames(res[[k]])<-NULL
			#An allele "dropps" if rfu <50
			res[[k]]$Heights<-replace(res[[k]]$Heights,which(is.na(res[[k]]$Heights)),0)
			#now that the matrix is built, we can easily modify it
			#vec  matches the observed and expected alles in the matrix, non observed (lost) alleles are coded NA
			vec<-rep(NA,length(res[[k]]$Expected))
			names(vec)<-res[[k]]$Expected
			locV<-as.character(intersect(res[[k]]$Expected,res[[k]]$Observed))
			vec[locV]<-locV
			#the same thing is done with peak heights
			He<-rep(0,length(res[[k]]$Expected))
			names(He)<-res[[k]]$Observed
			hh<-res[[k]]$Heights
			names(hh)<-res[[k]]$Observed
			He[which(res[[k]]$Expected %in% locV)]<-hh[locV]
			res[[k]]$Heights<-He
			res[[k]]$Observed<-vec
			Dtmp<-rep(0,length(res[[k]]$Expected))
			Dtmp2<-rep(0,length(res[[k]]$Expected))
			#the dropout variable is constructed by recording the alleles not observed
			for(h in 1:length(res[[k]]$Expected))
			{
				if(is.na(res[[k]]$Observed[h]) || res[[k]]$Heights[h]<s)  {Dtmp[h]<-1} #allele has dropped out
				else{ Dtmp[h]<-0} #no drop out
			}
			res[[k]]$D<-Dtmp
		}
		res

	}
	#if the data is a mixture, then the procedure must be modified accordingly
	else
	{
		for(k in 1:length(loci)) 
		{	
			#for each locus, store the allele counts for each person
			mat<-as.matrix(allele.count(genotype(as.character(as.matrix(geno[c(x,y),k])))))
			if(loci[k]!="AMEL" & loci[k]!="Amel")
			{
				if(ncol(mat)!=1)
				{
					tmp<-sort(as.numeric(colnames(mat)))
					tmp<-as.character(tmp)
					mat3<-t(as.matrix(mat[,tmp]))
				}
				if(ncol(mat)==1)
				{
					mat3<-t(mat)
				}
			}
			if(loci[k]=="AMEL" || loci[k]=="Amel")
			{
				if(ncol(mat)==1)
				{
					mat3<-t(mat)
				}
				if(ncol(mat)!=1)
				{
					tmp<-sort(colnames(mat))
					mat3<-t(as.matrix(mat[,tmp]))
				}
			}
			
			#temporary matrix storing the results for allele counts expected and observed alleles and the peak heights
			restmp<-cbind.data.frame(mat3,t(tabcsv[tabcsv$Marker==loci[k],c("Allele1","Allele2","Allele3","Allele4")][1:ncol(mat)]),t(tabcsv[tabcsv$Marker==loci[k],c("Height1","Height2","Height3","Height4")][1:ncol(mat)]))
			res[[k]]<-cbind.data.frame(rownames(restmp),restmp)
			colnames(res[[k]])<-c("Expected",paste("c",x,sep=""), paste("c",y,sep=""),'Observed','Heights')
			rownames(res[[k]])<-NULL
			res[[k]]$Heights<-replace(res[[k]]$Heights,which(is.na(res[[k]]$Heights)),0)
			#now that the matrix is built, we can easily modify it
			#vec  matches the observed and expected alleles in the matrix, non observed (lost) alleles are coded NA
			vec<-rep(NA,length(res[[k]]$Expected))
			names(vec)<-res[[k]]$Expected
			locV<-as.character(intersect(res[[k]]$Expected,res[[k]]$Observed))
			vec[locV]<-locV
			#the same thing is done with peak heights
			He<-rep(0,length(res[[k]]$Expected))
			names(He)<-res[[k]]$Observed
			hh<-res[[k]]$Heights
			names(hh)<-res[[k]]$Observed
			He[which(res[[k]]$Expected %in% locV)]<-hh[locV]
			res[[k]]$Heights<-He
			res[[k]]$Observed<-vec
			Dtmp<-rep(0,length(res[[k]]$Expected))
			Dtmp2<-rep(0,length(res[[k]]$Expected))
			#the dropout variable is constructed by recording the alleles not observed 
			for(h in 1:length(res[[k]]$Expected))
			{
				if(is.na(res[[k]]$Observed[h]) || res[[k]]$Heights[h]<s)  {Dtmp[h]<-1} #allele has dropped out
				else{ Dtmp[h]<-0} #no drop out
			}
			res[[k]]$D<-Dtmp
		}
		res
	}

}