dropDB <-
function(dataBase,DNA,dropHetero=0.05,alpha=0.5,dropIn=0.05,rel=c(0,0)/4,maxUnknown=1,adj=0,fst=0){

# Copyright (C) 20010 David J Balding, University College London.

# Version 1.0	19/1/10
# Version 1.1	23/1/10 restructured code to make it more similar to LR2unk.R
# Version 1.1.1	26/3/10 fixed small bug reported by Kirk Lohmueller, affecting the assignment of pps in 3 places
# Version 1.1.2	24/5/10 changed way dropin is modelled: only makes a difference if cc is not very small

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details see <http://www.gnu.org/licenses/>.

# The file LR1unk.R gives R functions for evaluating numerator and denominator of the single-locus LR for the hypothesis that a given, profiled individual (ss) is the unknown contributor to the crime scene profile (csp), versus another unprofiled individual (who may be related to ss).
# We assume exactly one unknown contributor but there can be any number of known contributors (whose alleles are masking alleles).
# Stutter alleles can also be included in the set of masking alleles: this allows the possibility that an apparent stutter could mask an allele.
# Multiple replicates are accommodated; currently dropout and dropin rates must be the same for each replicate, but this could be altered.
# There is no upper limit on the number of dropout or dropin alleles; it would be unusual to have more than one or two dropins at any locus, or more than say three or four in the profile; otherwise the hypothesis of another unknown contributor may be more realistic than multiple dropins; in practice cc is small and so penalises multiple dropins.
# Alleles are passed as names (character strings) and converted into indicator vectors

# For further details see the paper:

# Balding DJ, Buckleton J, Interpreting low template DNA profiles, Forensic Science International: Genetics, 4: 1-10, 2009, doi: 10.1016/j.fsigen.2009.03.003

# Relatedness was not discussed in the FSIG paper, this is a later addition. 

p = function(n,sub,adj,fst){
	
# Function p returns allele proportions after sampling and Fst adjustments.  Inputs are a vector n of integer allele counts and a vector sub that identifies the alleles to be corrected upwards; remaining alleles are corrected downwards, and the output is a non-negative vector that sums to 1.
# Usually sub includes the alleles of ss, and so the Fst adjustment allows for shared ancestry of ss with an alternative possible culprit;  if length(sub)=1 then a double adjustment is made (ss homozygous).
# I recommend Fst should be at least 0.01, and may need to be as high as 0.05 in some populations (e.g. small, isolated subpopulations of the population from which the reference database has been drawn).

		if(length(sub)==1) n[sub] = n[sub]+2*adj   # sampling adjustment
		else for(i in 1:length(sub)) n[sub[i]] = n[sub[i]] + adj
		n = n/sum(n)*(1-fst)/(1+fst)  # Fst adjustment
		if(length(sub)==1) n[sub] = n[sub]+2*fst/(1+fst)
		else for(i in 1:length(sub)) n[sub[i]] = n[sub[i]]+fst/(1+fst)
		n
}

kdrop = function(a,d,k){ # returns dropout prob for k copies of an allele when the prob for a single dropout is d
	a^(k-1)*d^k # this is a generalisation to more than two alleles of the homozygote dropout probability suggested by Balding & Buckleton (2009); there is no consensus about the best way to model multi-allele dropout, and the function given here can be replaced with another function if preferred; this function performs poorly if D is very large (close to 1) and a < 1
}

LRnumer1 = function(ss,call,n,d=0.2,a=0.5,cc=0.02){
	
# output: numerator of LR.

# inputs:
# ss specifies the profile of the alleged contributor; homozygotes should be coded as a scalar, heterozygotes as a vector of length 2.
# call is a list of lists.  The length of the outer list is the number of replicates (no upper limit).  
#Each inner list is of length 2: [[1]] crime scene profile (csp) alleles; [[2]] masking alleles (msk); 
#these are alleles observed in the profile that can be attributed to a known contributor #
#or to stutter that could mask the allele of an unknown contributor.  The program assumes that csp and msk have no allele in common, but there is currently no check for this.
# d, a, and cc control the dropout rates for heterozygotes and homozygotes, and the drop-in rate.

	pp = p(n,ss,adj,fst); # turn allele counts into fractions with sampling and Fst adjustment 
	ccsum = (1-cc^(1+length(pp)))/(1-cc) # sum of dropin probabilities; number of dropin alleles has geometric distribution with mean cc, but truncated at length(pp) so cc is not exactly the mean - but very close if cc is small
	csp = matrix(0,length(call),length(pp)) # initiate matrix of indicators of csp alleles, 1 row per replicate, 1 column per allele for which a database count is recorded
	msk = matrix(0,length(call),length(pp)) # initiate matrix of indicators of masking alleles, 1 row per replicate, 1 column per allele
	for(z in 1:length(call)){ # z indexes replicates
		tmp = call[[z]][[1]] 
		if(length(tmp)>0) for(i in 1:length(tmp)) csp[z,] = csp[z,] + (names(pp) == tmp[i]) 
		tmp = call[[z]][[2]] 
		if(length(tmp)>0) for(i in 1:length(tmp)) msk[z,] = msk[z,] + (names(pp) == tmp[i])
	}
	if(length(ss)==1) ss = c(ss,ss)  # convert scalar homozygote into a vector length 2
	xss = (names(pp)==ss[1]) + (names(pp)==ss[2]) # encode the two alleles of ss as a vector
	lrn = 1
	for(z in 1:length(call)){ 
		for(u in 1:length(pp)) if((xss[u]>0) & (msk[z,u]==0)){
			if(csp[z,u]==0) lrn = lrn * kdrop(a,d,xss[u])  # contribution from dropouts
			else lrn = lrn * (1-kdrop(a,d,xss[u])) # contribution from non-dropouts
		}
		din = csp[z,]*(xss==0) # vector indicating dropin alleles
		pps = pp/(1-sum(pp*msk[z,])) # allele proportions conditional on non-masking
		lrn = lrn * cc^sum(din)/ccsum * prod(pps[(1:length(pp))*din]) # contribution from dropins (= 1/ccsum if no dropins) 
	}
	lrn
}

LRdenom1 = function(ss,call,rel=c(0,0),n,d=0.2,a=0.5,cc=0.02){
# evaluates denominator of LR for arbitrary relationship between ss and alternative possible contributor.  

# Inputs the same as for LRnumer.

# Similar to LRnumer except that all possible genotypes are considered, instead of the genotype of the alleged contributor ss and the term for each genotype is weighted by an estimate of its population proportion, p^2 or 2pipj (assuming HWE).  If there is relatedness, then we consider the possibility that an allele is shared ibd with ss, in which case we only need to account for the other allele.

# First consider case of 0 alleles ibd
	pp = p(n,ss,adj,fst);
	ccsum = (1-cc^(1+length(pp)))/(1-cc)
	csp = matrix(0,length(call),length(pp)); msk = matrix(0,length(call),length(pp))
	for(z in 1:length(call)){
		tmp = call[[z]][[1]]
		if(length(tmp)>0) for(i in 1:length(tmp)) csp[z,] = csp[z,] + (names(pp) == tmp[i]) 
		tmp = call[[z]][[2]]
		if(length(tmp)>0) for(i in 1:length(tmp)) msk[z,] = msk[z,] + (names(pp) == tmp[i])
	}
	if(length(ss)==1) ss = c(ss,ss)
	lrd = 0
	for(i in 1:length(pp)) for(j in i:length(pp)){ # i and j specify a possible genotype for the alternative contributor
		xij = rep(0,length(pp))
		xij[i] = 1; xij[j] = xij[j]+1; # xij represents genotype ij as a vector
		term = pp[i]*pp[j]; if(i != j) term = 2*term # estimate of population proportion of genotype ij
		for(z in 1:length(call)){
			for(u in 1:length(pp)) if((xij[u]>0) & (msk[z,u]==0)){
				if(csp[z,u]==0) term = term * kdrop(a,d,xij[u])
				else term = term * (1-kdrop(a,d,xij[u]))
			}
			din = csp[z,]*(xij==0) # vector indicating dropin alleles
			pps = pp/(1-sum(pp*msk[z,])) 
			term = term * cc^sum(din)/ccsum * prod(pps[(1:length(pp))*din])
		}
		lrd = lrd + term
	}
	lrd = lrd * (1-sum(rel))  # multiply by probability of 0 alleles ibd

# Now 1 allele ibd
	if(rel[1] > 0) for(i in 1:2) for(j in 1:length(pp)){
		xsj = names(pp)==ss[i]  # ith allele of ss is the ibd allele
		xsj[j] = xsj[j]+1 # allele j makes up the genotype of the alternative contributor (related to ss)
		term = pp[j]
		for(z in 1:length(call)){
			for(u in 1:length(pp)) if((xsj[u]>0) & (msk[z,u]==0)){
				if(csp[z,u]==0) term = term * kdrop(a,d,xsj[u])
				else term = term * (1-kdrop(a,d,xsj[u]))
			}
			din = csp[z,]*(xsj==0)
			pps = pp/(1-sum(pp*msk[z,])) 
			term = term * cc^sum(din)/ccsum * prod(pps[(1:length(pp))*din]) 
		}
		lrd = lrd + rel[1]/2 * term  # rel[1] is divided by 2 because of two possibilities for the ibd allele
	}
	as.numeric(lrd)
}


LRnumer2 = function(ss,call,n,d=0.2,a=0.5,cc=0.02){
	
# output: numerator of LR.  See LR1unk.R for explanation of inputs and for comments on code

	pp = p(n,ss,adj,fst); # turn allele counts into fractions with sampling and Fst adjustment 
	ccsum = (1-cc^(1+length(pp)))/(1-cc) # sum of dropin probabilities, no. dropins has geometric distribution with parameter cc (truncated at length(pp))
	csp = matrix(0,length(call),length(pp)) # matrix of indicators of csp alleles, 1 row per replicate
	msk = matrix(0,length(call),length(pp)) # matrix of indicators of masking alleles, 1 row per replicate
	for(z in 1:length(call)){ # z indexes replicates
		tmp = call[[z]][[1]] 
		if(length(tmp)>0) for(i in 1:length(tmp)) csp[z,] = csp[z,] + (names(pp) == tmp[i]) 
		tmp = call[[z]][[2]] 
		if(length(tmp)>0) for(i in 1:length(tmp)) msk[z,] = msk[z,] + (names(pp) == tmp[i])
	}
	if(length(ss)==1) ss = c(ss,ss)  # convert scalar homozygote into a vector length 2
	lrn = 0
	for(i in 1:length(pp)) for(j in i:length(pp)){ # i and j specify a possible genotype of the one unknown contributor (ss is assumed to be a contributor)
		xijss = (names(pp)==ss[1]) + (names(pp)==ss[2]) # encode the two alleles of ss as a vector
		xijss[i] = xijss[i]+1; xijss[j] = xijss[j]+1  # add one to xijss for each of i and j
		term = pp[i]*pp[j]; if(i != j) term = 2*term # estimate of population fraction for genotype ij
		for(z in 1:length(call)){
			for(u in 1:length(pp)) if((xijss[u]>0) & (msk[z,u]==0)){
				if(csp[z,u]==0) term = term * kdrop(a,d,xijss[u])
				else term = term * (1-kdrop(a,d,xijss[u]))
			}
			din = csp[z,]*(xijss==0)
			pps = pp/(1-sum(pp*msk[z,])) 
			term = term * cc^sum(din)/ccsum * prod(pps[(1:length(pp))*din])
		}
		lrn = lrn + term
	}
	lrn
}

LRdenom2 = function(ss,call,rel=c(0,0),n,d=0.2,a=0.5,cc=0.02){

# First consider case that both unknown contributors have 0 alleles ibd with ss
	pp = p(n,ss,adj,fst)
	ccsum = (1-cc^(1+length(pp)))/(1-cc)
	csp = matrix(0,length(call),length(pp)) # matrix of indicators of csp alleles, 1 row per replicate
	msk = matrix(0,length(call),length(pp)) # matrix of indicators of masking alleles, 1 row per replicate
	for(z in 1:length(call)){ # z indexes replicates
		tmp = call[[z]][[1]] 
		if(length(tmp)>0) for(i in 1:length(tmp)) csp[z,] = csp[z,] + (names(pp) == tmp[i]) 
		tmp = call[[z]][[2]] 
		if(length(tmp)>0) for(i in 1:length(tmp)) msk[z,] = msk[z,] + (names(pp) == tmp[i])
	}
	if(length(ss)==1) ss = c(ss,ss)  # convert scalar homozygote into a vector length 2
	lrd = 0
	for(i in 1:length(pp)) for(j in i:length(pp)) for(k in 1:length(pp)) for(m in k:length(pp)){ # ij and km specify genotypes for the two unrelated unknown contributors 
		xijkm = rep(0,length(pp))
		xijkm[i] = 1; xijkm[j] = xijkm[j]+1; xijkm[k] = xijkm[k]+1; xijkm[m] = xijkm[m]+1  # xijkm gives the counts of the alleles for both unknown contributors
		term = pp[i]*pp[j]*pp[k]*pp[m]; if(i != j) term = 2*term; if(k != m) term = 2*term # estimate of population fraction for genotypes ij and km
		for(z in 1:length(call)){
			for(u in 1:length(pp)) if((xijkm[u]>0) & (msk[z,u]==0)){
				if(csp[z,u]==0) term = term * kdrop(a,d,xijkm[u])
				else term = term * (1-kdrop(a,d,xijkm[u]))
			}
			din = csp[z,]*(xijkm==0) # vector indicating dropin alleles
			pps = pp/(1-sum(pp*msk[z,])) 
			term = term * cc^sum(din)/ccsum * prod(pps[(1:length(pp))*din])
		}
		lrd = lrd + term
	}
	lrd = lrd * (1-sum(rel))

# Now the first unknown contributor has 1 allele ibd with ss
	if(rel[1] > 0) for(i in 1:2) for(j in 1:length(pp)) for(k in 1:length(pp)) for(m in k:length(pp)){
		xsjkm = names(pp)==ss[i] # choose an allele of s to be the ibd allele
		xsjkm[j] = xsjkm[j]+1; xsjkm[k] = xsjkm[k]+1; xsjkm[m] = xsjkm[m]+1
		term = pp[j]*pp[k]*pp[m]; if(k != m) term = 2*term
		for(z in 1:length(call)){
			for(u in 1:length(pp)) if((xsjkm[u]>0) & (msk[z,u]==0)){
				if(csp[z,u]==0) term = term * kdrop(a,d,xsjkm[u])
				else term = term * (1-kdrop(a,d,xsjkm[u]))
			}
			din = csp[z,]*(xsjkm==0) # vector indicating dropin alleles
			pps = pp/(1-sum(pp*msk[z,])) 
			term = term * cc^sum(din)/ccsum * prod(pps[(1:length(pp))*din])
		}
		lrd = lrd + rel[1]/2 * term
	}
	lrd
}


###Balding's code is above

nMarkers=length(unique(DNA$Marker))
mixture=known=defendant=mask=csp=missing=vector("list",nMarkers)
d1=dim(DNA)[1]
for (i in 1:nMarkers){
  mixture[[i]]=DNA[4*(i-1)+1,-c(1:2)]
  mixture[[i]]=mixture[[i]][!is.na(mixture[[i]])]
  known[[i]]=DNA[4*(i-1)+2,3:4]
  defendant[[i]]=DNA[4*(i-1)+3,3:4]
  missing[[i]]=DNA[4*(i-1)+4,-c(1:2)]
  missing[[i]]=missing[[i]][!is.na(missing[[i]])]
  csp[[i]]=setdiff(mixture[[i]],known[[i]])
  mask[[i]]=intersect(known[[i]],mixture[[i]])
}

marker.names=unique(dataBase$marker)
LR=numerator=denominator=NULL
for (i in 1:nMarkers){
 marker=marker.names[i]
 db=dataBase[dataBase$marker==marker,]
 n=db$count
 names(n)=db$allele
 alleles.defendant=as.character(defendant[[i]])
 alleles.csp.mask=list(list(csp[[i]],mask[[i]],missing[i]))
 if (maxUnknown==1){
  nmr = LRnumer1(alleles.defendant,alleles.csp.mask,n,dropHetero,alpha,dropIn)
  den = LRdenom1(alleles.defendant,alleles.csp.mask,rel,n,dropHetero,alpha,dropIn)
  this.LR = nmr/(nmr*rel[2]+den)
 }
 if (maxUnknown==2){
  nmr = LRnumer2(alleles.defendant,alleles.csp.mask,n,dropHetero,alpha,dropIn)
  den = LRdenom2(alleles.defendant,alleles.csp.mask,rel,n,dropHetero,alpha,dropIn)
  this.LR = nmr/(nmr*rel[2]+den)
 }
  
numerator=c(numerator,nmr)
denominator=c(denominator,den)
LR=c(LR,this.LR)
}

res=data.frame(marker=c(as.character(marker.names),"product"),LR=c(round(LR,2),prod(LR)))
res=data.frame(res,num=c(numerator,prod(numerator)),den=c(denominator,prod(denominator)))
res
}

