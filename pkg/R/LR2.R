#likelihood ratios new version allowing empty samples
#probbility of the evidence under H, function LRmix

# first define the  function which allows computing the likelihood of a given hypothesis

likEvid<-function(Repliste,T,V,x,theta,prDHet,prDHom,prC,freq)
{	
	sortieR<-0


	appelC<-function(Repliste,T,V,x,theta,prDHet,prDHom,prC,freq,sortieR)
	{

		lenRepliste<-length(Repliste)
		lenT<-length(T)
		lenV<-length(V)
		lenHom<-length(prDHom)
		lenHet<-length(prDHet)
		allele<-names(freq)
		freq0<-as.numeric(freq)
		lenFreq<-length(freq)
		sortieR<-0
		.C('evidenceC', as.double(Repliste), as.integer(lenRepliste),as.double(T), as.integer(lenT), as.double(V), as.integer(lenV),
		as.integer(x), as.double(theta), as.double(prDHet),as.integer(lenHet),as.double(prDHom),as.integer(lenHom),as.double(prC),as.character(allele),as.double(freq0),as.integer(lenFreq), sortieR=as.double(sortieR),PACKAGE='forensim')
		
		
	}
	
	tmp=(appelC(Repliste,T,V,x,theta,prDHet,prDHom,prC,freq,sortieR))
	return(tmp$sortieR)
}


LR2<-function(Repliste,Tp,Td,Vp,Vd,xp,xd,theta,prDHet,prDHom,prC,freq){

num<-likEvid(Repliste,T=Tp,V=Vp,x=xp,theta,prDHet,prDHom,prC,freq)
deno<-likEvid(Repliste,T=Td,V=Vd,x=xd,theta,prDHet,prDHom,prC,freq)


list('num'=num,'deno'=deno,'LR'=num/deno)}
