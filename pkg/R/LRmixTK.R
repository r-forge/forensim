# Hinda, Hague, june 2011
#GUI for LRmix
LRmixTK <-function()
{
	data(strusa)
	data(ngm)
	data(sgmNorway)
	#------------- Begin> formatting functions
	# function to convert tables of data, as exported from Genemapper to
		ConvertSamp<-function(tab)
		{
			whichsamp<-unique(tab[,1])#tclavle(sampname)# a verifier avec Corina
			whichmark<-unique(tab$Marker)

			if(any('AMEL' %in% whichmark)) whichmark<-whichmark[-which(whichmark=='AMEL')]
			replist<-vector('list',length(whichmark))
			names(replist)<-whichmark
			for(i in 1:length(whichmark))
			{
				reptmp<-NULL
				for(k in whichsamp)
				{		
					allele<-tab[tab$Marker==whichmark[i] & tab$SampleName==k,-c(1,2)]
					# print(allele)
					allele2<-as.numeric(allele[which(!is.na(allele))])
					# print('ici')

					
					list0<-list(allele2)
					# print(list0)
					names(list0)<-k
					reptmp<-c(reptmp, list0)
				}
				

				# print(reptmp)
				replist[[i]]<-reptmp
		}	
			
			replist
		}

	ConvertSamp2 <-
	function(tab)
	{
	whichsamp<-unique(tab[,1])#tclavle(sampname)# a verifier avec Corina
	whichmark<-unique(tab$Marker)

	if(any('AMEL' %in% whichmark)) whichmark<-whichmark[-which(whichmark=='AMEL')]
	replist<-vector('list',length(whichmark))
	names(replist)<-whichmark
	for(i in 1:length(whichmark))
	{
		reptmp<-NULL
		for(k in whichsamp)
		{		
			allele<-tab[tab$Marker==whichmark[i] & tab$SampleName==k,-c(1,2)]
			# print(allele)
			allele2<-(allele[which(!is.na(allele))])
			# print(class(allele2))
			a<-allele2[,1]
			b<-allele2[,2]
			
			
			list0<-list(paste(a,b,sep='/'))
			# print(list0)
			names(list0)<-k
			reptmp<-c(reptmp, list0)
		}
		# print(reptmp)
		replist[[i]]<-reptmp
	}	
		replist
	}


		#-------------End> formatting functions

	#Analyse function 
	analyse <-function(loc,repl, file1,file2,file3,ext1,ext2,ext3,infoStain)
	{

		# print(file3)
		
		# verify that the user explored the profiles and selected the loci 
		if(setequal(tclvalue(loc),''))
		{
			stop(tkmessageBox(message="First select the loci",icon="error",type="ok"))
			
		}
		else
		{
		# print(tclvalue(loc))
		loc0<- strsplit(tclvalue(loc),' ')[[1]]

		}
		
		# check replicates

		if(setequal(tclvalue(repl),''))
		{
			stop(tkmessageBox(message="First select the replicates",icon="error",type="ok"))
			
		}
		else
		{
		repl0<- strsplit(tclvalue(repl),' ')[[1]]

		}
		#check wether the user uploaded the files for the sample and the refrenecs profiles
		veriFile<-function(filename,ext,error=TRUE,txt='')
		{
				
					if(tclvalue(filename)=='')
					{	
						if(error){					
						stop(tkmessageBox(message=paste("First load the",txt,sep="","profile"),icon="error",type="ok"))
						}
						else {return(0)}
						# else
							# tkmessageBox(message="First load the reference profile",icon="info",type="ok")
					}
					else
					{
						if(tclvalue(ext)=='txt')
						{
							tab<-read.table(tclvalue(filename),h=TRUE,as.is=TRUE,sep='\t',na.string='')
							tab<-tab[-which(tab$Marker=='AMEL'),]

							#strings as strings, avoid converting to factors
						}
						else
						{
							tab<-read.csv(tclvalue(filename),h=TRUE,as.is=TRUE,na.string='')#strings as strings, avoid converting to factors
							tab<-tab[-which(tab$Marker=='AMEL'),]

						}# rm(file
						return(tab)
					}
				
		}

		
		
		
		
		#suspect samplename to be inserted into the listbox of the contributors under Hp, or eventually the non contributors under Hd
		
			verifFormat<-function(tab)
			{
				if(infoStain[1]!='SampleName')
				{
					tkmessageBox(message="Format error, pleae check your file",icon="error",type="ok")
				}
				if(infoStain[2]!='Marker')
				{
				# print(infoStain[2])
					tkmessageBox(message="Format error, please check your file",icon="error",type="ok")
				}
			}
			
			
			# CSP
			# if(tclvalue(extens1)=='txt')
			# stainFile<-read.table(tclvalue(filePath),h=TRUE,as.is=TRUE,sep='\t',na.string='')#strings as strings, avoid converting to factors
			# else{
			# stainFile<-read.csv(tclvalue(filePath),h=TRUE,as.is=TRUE,na.string='')#strings as strings, avoid converting to factors
			# }
		#suspect
		if(!require(tcltk)) stop("package tcltk is required")
		if(!require(tcltk2)) stop("package tcltk2 is required")
		if(!require(tkrplot)) stop("package tkrplot is required")
		tclRequire("Tktable")
		font0 <- tkfont.create(family="courrier",size=35,weight="bold",slant="italic")
		font1<-tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
		font2<-tkfont.create(family="times",size=16,weight="bold",slant="italic")
		font3<-tkfont.create(family="times",size=12)#,slant="italic")
		font4<-tkfont.create(family="courrier",size=14)#,slant="italic")
		font5<-tkfont.create(family="courrier",size=13,weight="bold")#,slant="italic")
		font6<-tkfont.create(family="times",size=8,weight='bold')#tkframe entries labels
		
		main <- tktoplevel()
		main2 <- tkframe(main)
		# tkgrid(main2)
		frame.tit<-tkframe(main2,relief="groove", borderwidth=4)
		tkgrid(frame.tit)
		tkgrid(tklabel(frame.tit,text="LR calculation", font="courrier 22", 	foreground="darkblue"),sticky='n')
		tkwm.title(main,"Analyse the profiles")

	#readFile button
			# nameInput <- tclvalue(tkgetOpenFile(parent=msmain,initialdir=tclvalue(path),multiple="true",
		firstFrame<-tkframe(main, relief="groove", borderwidth=4)
		tkgrid(tklabel(firstFrame,text='Parameters',font='courrier 14', fg='blue'))
		
		ncFrame<-tkframe(firstFrame,relief="flat", borderwidth=4)
		secondFrame<-tkframe(main,relief="groove", borderwidth=4)

		#---first frame 
		probFrame <- tkframe(firstFrame, relief="flat", borderwidth=4)
		repFrame<- tkframe(firstFrame, relief="groove", borderwidth=4)
		#--- second frame
		# extraFrame<-tkframe(secondFrame,relief='groove', borderwidth=4)
		hypoFrame<-tkframe(secondFrame, relief="groove", borderwidth=1)
		
		#third frame for uploading allele frequencies
		thirdFrame<-tkframe(main,relief="groove", borderwidth=4)
		
		
		#bottom frame
		bottomFrame<-tkframe(main,relief='flat',borderwidth=4)
		
		
		#-------------------- number of contributors ----------------- #
		tkgrid(tklabel(ncFrame,text="   Unknown contributors     ",font='courrier 10 bold'),sticky="w")
		ncHd<-tclVar(1)
		ncHp<-tclVar(1)
		
		ncHd.entry<-tkentry(ncFrame,textvariable=ncHd,width=4,highlightthickness=1,relief="solid",justify="center")
		ncHp.entry<-tkentry(ncFrame,textvariable=ncHp,width=4,highlightthickness=1,relief="solid",justify="center")
		tkgrid(tklabel(ncFrame,text=" Under Hp   ",font='courrier 10'),ncHp.entry,sticky="w")
		tkgrid(tklabel(ncFrame,text=" Under Hd   ",font='courrier 10'),ncHd.entry,sticky="w")
		

		
		#-------------------Hypotheses
		
		hypoFrame1<-tkframe(secondFrame)
		titre<-tkframe(hypoFrame1)
		titre2<-tkframe(thirdFrame)

		tkgrid(titre,pady=10)
		tkgrid(tklabel(titre,text="Hypotheses",font='courrier 14', fg='blue'))
		
		tkgrid(titre2,pady=10)
		tkgrid(tklabel(titre2,text="Allele Frequencies",font='courrier 14', fg='blue'))
		freqFrame<-tkframe(thirdFrame,relief="groove", borderwidth=4)
		#--- what kit to use????-----------------
		# tkgrid(tkcheckbutton(inFrame3, text="Diploid", variable=dip,font=font6, 
		# command=function() if (!as.logical(tclObj(dip))) tkconfigure(dip.entry, state="normal") else tkconfigure(dip.entry, state="disabled") ))
		sgm.var<-tclVar(0)
		sgmNorway.var<-tclVar(0)

		ngm.var<-tclVar(0)
		
		tkgrid(tkcheckbutton(thirdFrame, text='SGM+ US Caucasian', variable=sgm.var,
		font='courrier 8'), sticky='w')
		tkgrid(tkcheckbutton(thirdFrame, text='SGM+ Norwegian', variable=sgmNorway.var,
		font='courrier 8'), sticky='w')
		tkgrid(tkcheckbutton(thirdFrame, text='NGM', variable=ngm.var,font='courrier 8'),sticky='w')

		#--- first frame
		hypoFrame1.tit<-tklabel(hypoFrame1,text='Contributors under Hp',font='courrier 10 bold')
		hypoFrame11<-tkframe(secondFrame,relief='groove')
		# hypoFrame11.tit<-tklabel(secondFrame,text='Non-contributors under Hp')


		hypoFrame2<-tkframe(secondFrame)
		hypoFrame2.tit<-tklabel(hypoFrame2,text='Contributors under Hd',font='courrier 10 bold')
		# hypoFrame22.tit<-tklabel(secondFrame,text='Non-contributors under Hd')

		#-------------------------PROSECUTION-------------------------------------------------------------#
		#suspects variables
			
		csp<-veriFile(file1,ext1,error=TRUE,'crime scene profile')
		# print(csp)
		suspect<-veriFile(file2,ext2,error=TRUE,'suspect profile')
		suspectID<-unique(suspect$SampleName)
		# print(suspect)
		# contributors under Hp

		#victim and profiles not necesarily loaded, depends on the case
		victim<-veriFile(file3,ext3,error=FALSE,'')
		if(is.data.frame(victim))
		{
			vicID<-unique(victim$SampleName)
			for(v in 1:length(vicID))
			{
				#choice of contributors under Hp/Hd: victim profiles (if applicable)
				assign(paste('vicHp',v,sep=''), tclVar(1))
				assign(paste('vicHd',v,sep=''), tclVar(1))
				#dorpout for victims, if applicable

			}
		}
		tkgrid(hypoFrame1.tit)
		for(i in 1:length(suspectID))
		{
			# print(length(suspectID))
			# Tk entry for suspect(s) under Hp
			tkgrid(tklabel(hypoFrame1, text=suspectID[i],font='courrier 8'),sticky='w')
			# Tk entry for suspect(s) under Hd
		}
		
		
		tkgrid(hypoFrame1,padx=10)
		tkgrid(hypoFrame2.tit,pady=2)
		# for(i in 1:length(suspectID))
		# {
			# tkgrid(tkcheckbutton(hypoFrame2, text=suspectID[i], variable=get(paste('suspectHd',i,sep='')),font='courrier 8'),sticky='w')
		# }
		tkgrid(hypoFrame2)
		
		if(is.data.frame(victim)){
		for(i in 1:length(vicID))
		{
				tkgrid(tkcheckbutton(hypoFrame1, text=vicID[i], variable=get(paste('vicHp',i,sep='')),font='courrier 8'),sticky='w')
		
				tkgrid(tkcheckbutton(hypoFrame2, text=vicID[i], variable=get(paste('vicHd',i,sep='')),font='courrier 8'),sticky='w')
						
		}
		}	
		
		
		
		
		#---------prob of dropout and dropoin
			prD<-tclVar(0.1)
			prC<-tclVar(0.01)
			theta<-tclVar(0)
			#distribution
			prD.entry<-tkentry(probFrame,textvariable=prD,width=4,highlightthickness=1,relief="solid",justify="center")
			prC.entry<-tkentry(probFrame,textvariable=prC,width=4,highlightthickness=1,relief="solid",justify="center")
			theta.entry<-tkentry(probFrame,textvariable=theta,width=4,highlightthickness=1,relief="solid",justify="center")
			# merde=tkbutton(main)
			tkgrid(tklabel(probFrame,text="   Pr(D), Pr(C), theta     ",font='courrier 10 bold'))#,sticky="w")
			tkgrid(tklabel(probFrame,text="   Probability of Dropout   Pr(D)     ",font='courrier 10'),prD.entry,sticky="w")
			tkgrid(tklabel(probFrame,text="   Probability of Contamination   Pr(C)     ",font='courrier 10'),prC.entry,sticky="w")
			tkgrid(tklabel(probFrame,text="   Theta Correction (Fst)   ",font='courrier 10'),theta.entry,sticky="w")
		

		
		#-------- ANALYSIS OF RESULTS ---------- #
		
		analyse.loc<-function()
		{
				
			#----------VICTIM PROFILES  get users choice: victim
			#default values, eventually modified by whta follows
			Td.vic<-Tp.vic<-0
			if(is.data.frame(victim))
			{
				selecVicHp<-selecVicHd<-rep(0,length(vicID))#selecExt<-rep(0,length(extraID))
				for(j in 1:length(vicID))
				{
					selecVicHp[j]<-as.numeric(tclvalue(get(paste('vicHp',j,sep=''))))
					selecVicHd[j]<-as.numeric(tclvalue(get(paste('vicHd',j,sep=''))))
				}
				if(length(which(selecVicHp==1))!=0)
				{
					Tp.vic<-victim[victim$SampleName %in% vicID[which(selecVicHp==1)],]
				}
				else
				{
					Tp.vic<-NULL
				}
						
				if(length(which(selecVicHd==1))!=0)
				{
					Td.vic<-victim[victim$SampleName %in% vicID[which(selecVicHd==1)],]
					
				}
				else
				{
					Td.vic<-NULL
				}
			}
			
			#------------EXTRA PROFILES
			Vd.ext<-Vp<-0#default values, modified if users upload extra profiles
			
			#------------ SUSPECT PROFILES
		
			
			Tp.sus<-Vd.sus<-suspect
			
			#-------------------------------------------
			# assign the contributors under Hp: Tp
			if(is.data.frame(Tp.vic)) 
			{
			
				Tp<-rbind.data.frame(Tp.sus, Tp.vic)
			}
			else{
			# print('ici')
			Tp<-Tp.sus}
			
			

			#-------------------------------------------
			
			
			#assing non-contributors under Hd: Vd
			listDF<-list(Vd.sus,Vd.ext)
			tmp<-NULL
			for(d in listDF)
			{
				if(is.data.frame(d)) tmp<-rbind.data.frame(tmp,d)
			}
				
			# if(is.null(tmp)){ Vd<-NULL}
			# else{ Vd<-tmp}
			# print(Vd)
			cspFinal<-ConvertSamp(csp[csp$Marker %in% loc0 & csp$SampleName %in% repl0, ])
			#contributors under Hp
			
			TpFinal<-sapply(ConvertSamp2(Tp[Tp$Marker %in% loc0, ]), unlist)
			if(is.character(TpFinal)) TpFinal <-matrix(TpFinal,ncol=length(loc0))
			
			if(is.data.frame(Vp))
			{
				VpFinal<-sapply(ConvertSamp2(Vp[Vp$Marker %in% loc0, ]), unlist)
				if(is.character(VpFinal)) VpFinal <-matrix(VpFinal,ncol=length(loc0))

			}
			else{VpFinal<-Vp}
			
			# if(is.data.frame(Vd))
			# {
			# VdFinal<-sapply(ConvertSamp2(Vd[Vd$Marker %in% loc0, ]), unlist)
			# if(is.character(VdFinal)) VdFinal <-matrix(VdFinal,ncol=length(loc0))

			# }
			# else{VdFinal<-Vd}
			
			if(is.data.frame(Td.vic))
			{
			TdFinal<-sapply(ConvertSamp2(Td.vic[Td.vic$Marker %in% loc0, ]), unlist)
			if(is.character(TdFinal)) TdFinal <-matrix(TdFinal,ncol=length(loc0))

			}
			else{TdFinal<-Td.vic}
			# print(class(TpFinal))
			# print(ConvertSamp2(TpFinal))
			lr0<-rep(0,length(loc0))
			#if there are known non-contributors under Hp
			n1<-as.numeric(tclvalue(sgm.var))
			n2<-as.numeric(tclvalue(ngm.var))
			n3<-as.numeric(tclvalue(sgmNorway.var))

			if(n1==0 & n2==0 & n3==0)
			{
				stop(tkmessageBox(message="please choose the allele frequencies",icon="error",
				type="ok"))
			}
			if(length(which(c(n1,n2,n3) %in% 1))>1)
			{
				stop(tkmessageBox(message="only one choice is possible for the  allele
				frequencies",icon="error",type="ok"))
			}
			#if the user made the right correction
			
			if(n1==1){data0<-strusa$tab$Cauc}

			if(n2==1){data0<-ngm$tab}
			if(n3==1){data0<-sgmNorway$tab}

			for(jj in 1:length(loc0))
			{
				mark0<-loc0[jj]
				rep0<-cspFinal[[jj]]
				# if(is.matrix(VpFinal)){ tmpVp<-VpFinal[,jj]}
				# else{ tmpVp<-NULL}
				
				if(is.matrix(TdFinal )){ tmpTd<-TdFinal[,jj]}
				else{ tmpTd<-NULL}
				
				# if(is.matrix(VdFinal)){ tmpVd<-VdFinal[,jj]}
				# else{ tmpVd<-NULL}
			
				#numerator Pr(E|Hp)
				lr0[jj]<-LR2(Repliste=rep0,Tp=TpFinal[,jj],Vp=NULL,xp=as.numeric(tclvalue(ncHp)),
				Td=tmpTd,Vd=TpFinal[,jj],xd=as.numeric(tclvalue(ncHd)), 
				theta=as.numeric(tclvalue(theta)),prDHet=as.numeric(tclvalue(prD)),
				prDHom=as.numeric(tclvalue(prD))*as.numeric(tclvalue(prD))*0.5,
				prC=as.numeric(tclvalue(prC)),freq=data0[[mark0]])$LR
				
				
			}
			
			#display the results
			res<-tktoplevel()
			tkwm.title(res,"LRmix: Results            ")

			f1<-tkframe(res)
			tkgrid(tklabel(f1,text="   Results     ",font='courrier 14', foreground="blue"),sticky="w")
			
			#
			array1 <- tclArray()
			array2 <- tclArray()

			myRarray<-c("LR per Locus",loc0,"LR",signif(lr0,4))
			myRarray2<-c("Overall LR",signif(prod(lr0),4))

			dim(myRarray)<-c(length(loc0)+1,2)
			dim(myRarray2)<-c(2,1)

			for(i in 0:(length(loc0))){
			array1[[i,0]] <- myRarray[i+1,1]
			array1[[i,1]] <- myRarray[i+1,2]

			}
			# print(myRarray)
			#table2
			array2[[0,0]] <- myRarray2[1,1]
			array2[[1,0]] <- myRarray2[2,1]

			#if no known non-contributors under Hp     
			# print('length');print((loc0))
			table1 <- tkwidget(f1,"table",variable=array1,rows=(length(loc0))+1 ,cols=2,titlerows=1,titlecols=0, colwidth=25)
			
			table2 <- tkwidget(f1,"table",variable=array2,rows=2,cols=1,titlerows=1,titlecols=0, colwidth=25)
			
			
			#xscr <-tkscrollbar(tt,orient="horizontal", command=function(...)tkxview(table1,...))
			# yscr <- tkscrollbar(f1,repeatinterval=5, command=function(...)tkyview(table1,...))
			# tkgrid(table1, table2,yscr,sticky="new", padx=20,pady=18)
			tkgrid(table1, table2,sticky="new", padx=20,pady=18)
			# tkgrid.configure(yscr,sticky="nsw")
			#tkconfigure(table1,variable=array1,background="white",selectmode="extended")
			tkconfigure(table1,variable=array1,background="white",selectmode="extended",rowseparator="\"\n\"",colseparator="\"\t\"")
			tkconfigure(table2,variable=array2,background="white",selectmode="extended",rowseparator="\"\n\"",colseparator="\"\t\"")
			
			# button to display the plot of the LR vs the PrD
			plot.but<-tkbutton(f1, text="Plot LR vs PrD",fg="blue", font="courrier 10",command=function() Dplot())#,command=function() openFile())
			Dplot<-function()
			{
				LRres<-vector('list', length(loc0))
				vecD<-seq(0.00001,0.99,length=10)
				for(jj in 1:length(loc0))
				{
					tmp1<-rep(0,length(vecD))
					mark0<-loc0[jj]
					rep0<-cspFinal[[jj]]
					if(is.matrix(TdFinal )){ tmpTd<-TdFinal[,jj]}
					else{ tmpTd<-NULL}
				
					for(k in 1:length(vecD))
					{				
					#numerator Pr(E|Hp)
						tmp1[k]<-LR2(Repliste=rep0,Tp=TpFinal[,jj],Vp=NULL,
						xp=as.numeric(tclvalue(ncHp)),Td=tmpTd,Vd=TpFinal[,jj],
						xd=as.numeric(tclvalue(ncHd)),theta=as.numeric(tclvalue(theta)),
						prDHet=vecD[k],prDHom=vecD[k]*vecD[k]*0.5,prC=as.numeric(tclvalue(prC)),
						freq=data0[[mark0]])$LR
						
					}	
					#LR for locus jj
					LRres[[jj]]<-tmp1
				}
				tmp<-NULL
				for(i in LRres){
				tmp<-cbind(tmp,i)}

				tmp<-apply(tmp,1,prod)
				# print(tmp)
				if('Inf' %in% range(tmp) | '-Inf' %in% range(tmp)| NaN %in% range(tmp) ){
					
					stop(tkmessageBox(message="infinite LR values, please change the model parameters",icon="error",type="ok"))
				}
				
				
				Myhscale <- 1   # Horizontal scaling
				Myvscale <- 1
				dd <- tktoplevel()
				# tkconfigure(dd,cursor="watch")

				tkwm.title(dd,"LR plot")
				frameC<-tkframe(dd)
				# print('ici')
				
				# print(LRres)
				Dplot.loc<-function()
				{
					params <- par(bg="white")
					plot(vecD,log(tmp,10),ylab='log10 LR',xlab='Probability of Dropout',cex.lab=1.3,xlim=c(0,1),pch=19,ylim=range(log(tmp,10),finite=TRUE))
					lines(vecD, log(tmp,10),lty=3,col='gray')
					title('LR vs. probability of dropout', cex=1.3)
					par(params)
				}
				img <- tkrplot(frameC,fun= Dplot.loc,hscale=Myhscale,vscale=Myvscale)
				CopyToClip <- function()
				{
					tkrreplot(img)
				}
				copy.but <- tkbutton(frameC,text="Copy to Clipboard",font="courrier 10",fg="darkblue",command=CopyToClip)
				tkgrid(img)
				tkgrid(copy.but)
				tkpack(frameC)
				
			}
			
			
			
		# grid frame f1 contains the tables
		tkgrid(plot.but,pady=25,columnspan=10)
		tkgrid(f1)
		}

		
		
			
		go.butt<-tkbutton(bottomFrame,text='OK!', font='courrier 13',fg='blue', command= analyse.loc)#
		
		# 
		# print(head(Td.sus))

		#----- grid and pack, analysis frame
			tkgrid(hypoFrame1, pady=12,padx=10)

		
		tkgrid(ncFrame,pady=12)#,pady=10, padx=10)
		tkgrid(probFrame,pady=12)
		tkgrid(secondFrame,firstFrame, thirdFrame,pady=10,padx=10)#, pady=10, padx=10)

		
		# tkgrid(listHd, pady=15 )
		# tkgrid(listHd2,pady=10)
		# tkgrid(freqFrame)
		tkgrid(bottomFrame,pady=10,columnspan=45)
		tkgrid(go.butt,columnspan=45)

		

	}


	
	
	
	
	#main
	
	
	if(!require(forensim)) stop("package forensim is required")
	if(!require(gdata)) stop("package gdata is required")
	if(!require(gtools)) stop("package gtools is required")
	if(!require(MASS)) stop("package MASS is required")
	if(!require(mvtnorm)) stop("package mvtnorm is required")
	if(!require(genetics)) stop("package genetics is required")
	
	if(!require(tcltk)) stop("package tcltk is required")
	if(!require(tcltk2)) stop("package tcltk2 is required")
	if(!require(tkrplot)) stop("package tkrplot is required")
	tclRequire("Tktable")
	font0 <- tkfont.create(family="courrier",size=35,weight="bold",slant="italic")
	font1<-tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
	font2<-tkfont.create(family="times",size=16,weight="bold",slant="italic")
	font3<-tkfont.create(family="times",size=12)#,slant="italic")
	font4<-tkfont.create(family="courrier",size=14)#,slant="italic")
	font5<-tkfont.create(family="courrier",size=13,weight="bold")#,slant="italic")
	font6<-tkfont.create(family="times",size=8,weight='bold')#tkframe entries labels
	font7<-tkfont.create(family="courrier",size=10,weight='bold')#tkframe entries labels

	tf <- tktoplevel()
	tkwm.title(tf,"LRmix: Likelihood Ratio Calculator")
	# icn <- tkimage.create("photo", file=system.file("files/test.GIF", package = "forensim"))#"test.GIF")
	#TclTklabel <- tklabel(frame1, image=icn, background="white")
	done <- tclVar(0)
	filePath<-tclVar('')

	# filePath4<-tclVar('')
	extens1<-tclVar('')
	
	locus<-tclVar('')
	#replicates
	repl<-tclVar('')
	#data to use in the calculations
	# CSP<-tclVar()
	# suspect<-tclVar(0)
	# victim<-tclVar(0)
	# elimination<-tclVar(0)
	
	openFile<-function(file0,caselist,top,ext)
	{		
		fileName<-tclvalue(tkgetOpenFile(parent=top,initialdir=tclvalue(file0),multiple="true",
		filetypes="{{CSV Files} {.csv .txt}} {{Tab-delimited Files} {.tab}}")) #tclvalue(tkgetOpenFile())
		if (!nchar(fileName))
		{
			tkmessageBox(message="No file was selected!")
		}
		else
		{
			tmp<-sub('\\}',fileName,replacement='')
			tmp2<-sub('\\{',tmp,replacement='')
			tclvalue(file0)<-tmp2
			foo3<-strsplit(tmp2,'/')[[1]]
			tclvalue(ext)<-strsplit(foo3[length(foo3)],'\\.')[[1]][2]
			# print('here');print(foo3foo3[length(foo3)])
			tkinsert(caselist,0,paste(foo3[length(foo3)],sep=":"))
		}
	}
	
#---------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
'sampleProf'<-function()
{
	tprof <- tktoplevel()
	tkwm.title(tprof,"LRmix: import DNA sample profiles")
	filesFrame<-tkframe(tprof)
	tkgrid(tklabel(filesFrame, text="   DNA samples   ",font=font4,foreground="blue"), columnspan=9)
	bottomFrame<-tkframe(tprof)
	#white box will contain the lsist of the available files
	caselist <- tklistbox(filesFrame,height=10,selectmode="extended",background="white",width=25)
   #------------- Display profiles from crime stain, so that the user can choose the loci ---------#
   displayprofile<-function()
   {	
		prof<- tktoplevel()
		tkwm.title(prof,"Sample profile")
		done <- tclVar(0)
		f1 <- tkframe(prof, relief="groove", borderwidth=4)#,bg='white')
		repFrame1 <-tkframe(f1, relief="groove", borderwidth=4)#,bg='white')
		repFrame2 <-tkframe(f1, relief="groove", borderwidth=4)#,bg='white')
		repFrame3 <-tkframe(f1, relief="groove", borderwidth=4)#,bg='white')
		repFrame4 <-tkframe(f1, relief="groove", borderwidth=4)#,bg='white')

		f3 <-tkframe(f1, relief="groove", borderwidth=4)#,bg='white')
		f4 <-tkframe(prof, relief="flat", borderwidth=4)#,bg='white')
		tkgrid(tklabel(f1, text="      Loci selection        ",font=font4,foreground="blue"), columnspan=9)
		#a max of four replicates
		

		if(tclvalue(extens1)=='txt')
		stainFile<-read.table(tclvalue(filePath),h=TRUE,as.is=TRUE,sep='\t',na.string='')#strings as strings, avoid converting to factors

		else{
		stainFile<-read.csv(tclvalue(filePath),h=TRUE,as.is=TRUE,na.string='')#strings as strings, avoid converting to factors
		}
		
		#remove the Amel Marker
		stainFile<-stainFile[-which(stainFile$Marker=='AMEL'),]

		# print(tclvalue(CSP))
		#handling replicates
		# print(tclvalue(extens1))
		infoStain<-names(stainFile)
		# print(infoStain)
		# sample Infor;ation: replicates
		##### General verifications
		if(infoStain[1]!='SampleName')
		{
			tkmessageBox(message="Format error, pleae refer to the reference file",icon="error",type="ok")
		}
		if(infoStain[2]!='Marker')
		{
			# print(infoStain[2])
			tkmessageBox(message="Format error, pleae refer to the reference file",icon="error",type="ok")
		}
		
		
		# number of replicates
		sampname<-(unique(stainFile[,1]))
		
		if(length(sampname)>4)
		{
			tkmessageBox(message=paste("There are",length(sampname),"replicates, only the four first will be used"),icon="warning",	type="ok")
		}
		locusStain<-as.character(unique(stainFile$Marker))
		# print(locusStain)
		#okprof function selects the loci: these sould be selected only for the samples
		#and not the reference profiles
		
		okprof<-function()
		{
			selecLoci<-rep(0,length(locusStain))
			selecRep<-rep(0,length(sampname))
			for(k in 1:length(locusStain))
			{
				selecLoci[k]<-as.numeric(tclvalue(get(paste('loc',k,sep=''))))
				
			}
			
			for(k in 1:length(sampname))
			{
				selecRep[k]<-as.numeric(tclvalue(get(paste('tclRep',k,sep=''))))
				
			}
			#to get user's choice of replicates

			if(length(which(selecRep==1))==0)
			{
				# print(length(which(selecRep==1))==0)
				(tkmessageBox(message="At least one replicate must be selected",icon="error",type="ok"))
			}
			else {
			
				tclvalue(repl)<-sampname[which(selecRep==1)]
			}
			
			# get user's choice of the loci
			if(length(which(selecLoci==1))==0)
			{
				(tkmessageBox(message="At least one locus must be selected",icon="error",type="ok"))
			
			}
			
			else
			{
				tclvalue(locus)<-locusStain[which(selecLoci==1)]
			}

			# print(tclvalue(CSP))

			tkdestroy(prof)
			tkdestroy(tprof)
			# return(locusStain2)
		}
		
		#create tclVar locus
		for(m in 1:length(locusStain)){
		assign(paste('loc',m,sep=''), tclVar(1))}
		
		
		
		
		j=1
		while(j<(length(sampname)+1))#replicates
		{
			for(i in 1:length(locusStain))
			{
				tmp0<-stainFile[stainFile$SampleName==sampname[j],]
				# print(sampname[j])
				tmp1<-tmp0[tmp0$Marker==locusStain[i],-c(1,2)]
				tmpProf<-unlist(tmp1[which(!is.na(tmp1))])
				# print(tmpProf)
				names(tmpProf)<-NULL
				# print(tmpProf)
				tkgrid(tkcheckbutton(get(paste('repFrame',j,sep='')), text=locusStain[i], variable=get(paste('loc',i,sep='')),font='courrier 8'),sticky='w',columnspan=9)
				
				if(length(tmpProf)==0)#null profile
				{
					tkgrid(tklabel(get(paste('repFrame',j,sep='')), text=0),sticky='we', columnspan=9)
				}
				else
				{
					tkgrid(tklabel(get(paste('repFrame',j,sep='')), text=tmpProf),sticky='we', columnspan=9)
				}
				
			}
			j=j+1
		}

		
		tkpack(f1)
		# for(j in 1:length(sampname))
		# {
			# replicates[j]<-get(paste(rep,j,sep=''))
		# }	
		
		#--------------Add replicates checkbuttons
		#--------------Add replicates checkbuttons
		# first create the tcl variables
		for(m in 1:length(sampname)){
		assign(paste('tclRep',m,sep=''), tclVar(1))}
		for(i in 1:length(sampname))
		{
			tkgrid(tkcheckbutton(get(paste('repFrame',i,sep='')), text=sampname[i],fg='blue', variable=get(paste('tclRep',i,sep='')),font='courrier 14'),sticky='w',columnspan=9)
			# tkgrid(tkcheckbutton(get(paste('repFrame',i,sep='')), text=paste('Replicate',i,sep=''),fg='blue', variable=get(paste('tclRep',i,sep='')),font='courrier 14'),sticky='w',columnspan=9)
		}
		
		#second create the checkbuttons
		
		if(length(sampname)>=4)
		{
			tkgrid(repFrame1,repFrame2,repFrame3,repFrame4,padx=12,pady=8 )
			
			#,repFrame3,repFrame4,padx=12,pady=8 )

			# tkgrid(tklabel(rep1, text=paste("Replicate 1",sampname[1],sep=':'),font=font4,foreground="blue"), columnspan=9)
			# tkgrid(tklabel(rep2, text=paste("Replicate 2",sampname[2],sep=':'),font=font4,foreground="blue"), columnspan=9)
			# tkgrid(tklabel(rep3, text=paste("Replicate 3",sampname[3],sep=':'),font=font4,foreground="blue"), columnspan=9)
			# tkgrid(tklabel(rep4, text=paste("Replicate 4",sampname[4],sep=':'),font=font4,foreground="blue"), columnspan=9)
		}
		else
		{
		
			if(length(sampname)==1)
			{
				# tkgrid(tklabel(rep1, text=paste("Replicate 1",sampname,sep=':'),font=font4,foreground="blue"), columnspan=9)
				tkgrid(repFrame1,padx=12,pady=8 )
				
			}
			
			if(length(sampname)==2)
			{
				tkgrid(repFrame1,repFrame2,padx=12,pady=8 )
				# tkgrid(tklabel(rep1, text=paste("Replicate 1",sampname[1],sep=':'),font=font4,foreground="blue"), columnspan=9)
				# tkgrid(tklabel(rep2, text=paste("Replicate 2",sampname[2],sep=':'),font=font4,foreground="blue"), columnspan=9)
			
			}
			
			if(length(sampname)==3)
			{
			
				tkgrid(repFrame1,repFrame2,repFrame3,repFrame4,padx=12,pady=8 )
				# tkgrid(tklabel(rep1, text=paste("Replicate 1",sampname[1],sep=':'),font=font4,foreground="blue"), columnspan=9)
				# tkgrid(tklabel(rep2, text=paste("Replicate 2",sampname[2],sep=':'),font=font4,foreground="blue"), columnspan=9)
				# tkgrid(tklabel(rep3, text=paste("Replicate 1",sampname[3],sep=':'),font=font4,foreground="blue"), columnspan=9)
			
			}



		}
		
		okprof.but<-tkbutton(f4, text="OK!",fg="darkblue", font="courrier 12",command=function() okprof())#,command=function() openFile())
		# print(locus)
		tkgrid(okprof.but,pady=2)
		tkpack(f4)
   }

   
readFile.but<-tkbutton(bottomFrame, text="Import datafile",fg="darkblue", font="courrier 12",command=function() openFile(filePath,caselist,tf,extens1))
profile.but<-tkbutton(bottomFrame, text="Display profile",fg="darkblue", font="courrier 12",command=function() displayprofile())
tkgrid(filesFrame) 
tkgrid(caselist,padx=10,pady=10)
tkgrid(bottomFrame)
tkpack(readFile.but,  profile.but,side='left')
}		



#--------------------------------------------------------------------------------#
#--------- Reference profiles
#----------------------------------------------------------------------------------
filePath3<-tclVar(''); 	extens3<-tclVar('')
filePath2<-tclVar(''); 	extens2<-tclVar('')
'refProf'<-function()
{
	tref <- tktoplevel()
	tkwm.title(tref,"Reference profiles")
	#--------------------- refrence profiles---------------------#		
	refFrame<-tkframe(tref)
	buffer1<-tkframe(refFrame,relief='groove')
	buffer1.tit<-tkframe(refFrame,relief='groove')
	buffer2.tit<-tkframe(refFrame)
	buffer3.tit<-tkframe(refFrame)

	buffer2<-tkframe(refFrame)
	buffer3<-tkframe(refFrame,relief='groove')

	# tkgrid(tklabel(buffer1, text="   Suspect(s)   ",font=font4,foreground="blue"), columnspan=9)

	#------ suspect 
	suspectFrame<-tkframe(buffer1, relief='groove')

	#------- victim
	victimFrame<-tkframe(buffer2,relief='groove')
	elimFrame<-tkframe(refFrame)
	bottomFrame<-tkframe(refFrame)
	#-------- suspects frame scrollabr 
	scr2 <- tkscrollbar(buffer1, repeatinterval=10)#, command=function(...)tkyview(caselist,...))
	#white box will contain the lsist of the available files
	caselist2 <- tklistbox(buffer1,height=5,selectmode="extended",background="white",width=20)

	# tkgrid(refFrame)
	
	#-------victim(s) frame scrollabr 
	scr3 <- tkscrollbar(buffer2, repeatinterval=10)#, command=function(...)tkyview(caselist,...))
	#white box will contain the lsist of the available files
	caselist3 <- tklistbox(buffer2,height=5,selectmode="extended",background="white",width=20)
	
	
	#-------extra profiles frame scrollabr 
	# scr4 <- tkscrollbar(buffer3, repeatinterval=10)#, command=function(...)tkyview(caselist,...))
	#white box will contain the lsist of the available files
	# caselist4 <- tklistbox(buffer3,height=5,selectmode="extended",background="white",width=20)
	###---------Buttons: Open File and check?
	#openfile is called when button Openfile is pressed: values of filepath2,caselist2, etc are updated 
	ref.but<-tkbutton(suspectFrame, text="Open File",fg="darkblue", font="courrier 10", command= function() 
	openFile(filePath2,caselist2,tref,extens2))
	# check.but<-tkbutton(suspectFrame, text="Check?",fg="darkblue", font="courrier 10")#,command=function() openFile())
	# a lot of TCLvariables in order to 
	
	ref.but2<-tkbutton(victimFrame, text="Open File",fg="darkblue", font="courrier 10", command= function() 
	openFile(filePath3,caselist3,tref,extens3))
	# check.but2<-tkbutton(victimFrame, text="Check?",fg="darkblue", font="courrier 10")#,command=function() openFile())
	# a lot of TCLvariables in order to 
	
	# ref.but3<-tkbutton(elimFrame, text="Open File",fg="darkblue", font="courrier 10", command= function() 
	# openFile(filePath4,caselist4,tref,extens4))
	# check.but2<-tkbutton(victimFrame, text="Check?",fg="darkblue", font="courrier 10")#,command=function(
	
	
	ok.but<-tkbutton(bottomFrame,text=' OK ', font=font7,fg='blue',command=function()tkdestroy(tref))
	
	# checkInput<-function(){}
	#suspects
	tkgrid(caselist2)
	tkgrid(suspectFrame,padx=10,pady=10)
	tkgrid(buffer1.tit)
	tkgrid(tklabel(buffer1.tit, text="   Suspect(s)   ",font=font4,foreground="blue"), columnspan=9)
	tkgrid(buffer1,pady=10,padx=10)
	tkgrid(ref.but)#, check.but)
	# victims
	tkgrid(tklabel(buffer2.tit, text="   Victim(s)   ",font=font4,foreground="blue"),pady=3)
	tkgrid(caselist3)
	tkgrid(buffer2.tit)
	tkgrid(buffer2,pady=10,padx=10)
	tkgrid(victimFrame,pady=10,padx=10)
	# tkgrid(tklabel(buffer2, text="   Victim(s): known contributors   ",font=font4,foreground="blue"), columnspan=9)
	tkgrid(ref.but2)#, check.but2)
	tkgrid(ok.but,pady=2)
	tkgrid(bottomFrame)
	tkgrid(refFrame)	

}
	
	

	
	##############--------Main frame-----------------------------------------------
	
			
	main <- tkframe(tf, relief="groove", borderwidth=2)
	frame0<-tkframe(main, relief="groove", borderwidth=2)
	frame1<-tkframe(main)

	# frameButt <- tkframe(tf, relief="groove", borderwidth=4)
	#icn <- tkimage.create("photo", file=system.file("files/test.GIF", package = 		"forensim"))#"test.GIF")
	#TclTklabel <- tklabel(frame1, image=icn, background="white")
	labh <- tklabel(frame0,bitmap='questhead')#, image=icn)
	#labh <- tklabel(frame1)
	tkbind(labh, "<Button-1>", function() 'hh')
	# tkgrid(labh)
	
	
	#--------------
	# labh <- tklabel(tf)
	#labh <- tklabel(frame1)
	#tkbind(labh, "<Button-1>", function() 'hh')
	tkgrid(tklabel(frame0,text="   Evaluation of Likelihood Ratios  ", font="times 20", foreground="darkblue"),labh)
	
	#---------
	
	samp.butt<-	tkbutton(frame1, text=" Load Sample Profiles ",fg="darkblue", font="courrier 12", command=sampleProf)
	ref.butt<-tkbutton(frame1, text=" Load Reference Profiles ",fg="darkblue", font="courrier 12",command=refProf)
	# call analyse() with arguments locus and filepath, that are modified by 
	analyse.butt<-tkbutton(frame1, text=" Analyse ",fg="darkblue", font="courrier 12",command=function() analyse(locus,repl,filePath,filePath2,filePath3, extens1, extens2, extens3,infoStain))
	tkgrid(samp.butt, padx=10,pady=10)
	tkgrid(ref.butt,padx=10,pady=10)
	tkgrid(analyse.butt,padx=10,pady=10)
	tkgrid(frame0)
	tkgrid(frame1)
	tkgrid(main)

}

