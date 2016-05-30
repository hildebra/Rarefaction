#' Creates a rarefaction curve. Depending on the method, in some cases only richness histogram is returned
#' @param samples: number of samples at each rarefaction depth, used to average richness estimator
#' @param method: "S_obs", "rarefy", "shannon", "simpson","invsimpson","ACE","ICE","chao1","chao2"
#' ACE = abundance data, ICE = presence/absense
#' @param plotType box,single,vio
#' @param maxEstPlot box,density
#' @param recalcRare (T) forces to recalc the rarefied matrix, ignoring if the matrix was pre-calculated
Interface.AlphaDiversity = function(Mat,Opt=extractOptions(Mat),cls=NULL,
		MetaInfo=extractMetaInfo(Mat),h=NULL, samples = 1,method=.DDT.RichnessMethods,
		rarepoints = NULL, silent = FALSE, givenRichness = NULL,whatIs = Opt$what,
		plotType=c("box","single","vio"),subFeatures=NULL,subsection=TRUE,
		title="Taxa Richness estimation",DOtest=T,Cl=extractParClus(Mat),
		maxEstPlot= c("box","density"),subset=NULL,smplcol=NULL,poolSamples=NULL,
		recalcRare=F,extendRseries=T,removeSubRareSmpls=FALSE){
	#browser()
	myLibs("Biobase")
	cat("Starting richness / diversity analysis\n")
	ret 		= list();
	colsel      = NULL;
	class(ret) 	= "DDT.richness"
	plotType 	= match.arg(plotType)
	maxEstPlot 	= match.arg(maxEstPlot)
	method 		= match.arg(method)

	if (plotType%out%c("single")){extendRseries=FALSE}
	PHtree 		= extractPhyloTree(Mat,MetaInfo)
	if (is.null(PHtree) && method %in% c("FaithPD","MNTD","MPD")){
		writeInfo(paste("Error: no phylogenetic tree provided, can not caluclate",method),Opt,"warning")
		return(NULL);
	}
	#myLibs("gplots");
	#myLibs("vioplot");
	#preparation
	#if (!silent){writeInfo("",Opt,"mainFile",write=subsection);}
	writeInfo(title,Opt,"section",write=!silent);
	blck 		= attr(cls,"block")

	if(is.null(subFeatures)){subFeatures=c("all")
	}else if (!any(subFeatures=="all")){subFeatures = c("all",subFeatures)}
	if (is.null(subset)&&!is.null(cls) && !is.null(names(cls))){
		cls 	= cls[!is.na(cls)]
		subset 	= names(cls)
	}

	if(!is.null(givenRichness)){
		if (!is.null(subset)){
				mm 				= subsetSel(names(givenRichness),subset)
				givenRichness 	= givenRichness[mm]
		}
		nar=is.na(givenRichness)
		if (any(nar)){
			givenRichness 	= givenRichness[!nar]
			subset			= names(givenRichness)
			cls 			= cls[!nar];
		}
		subFeatures		= c("all")
		whatIs 			= ""
	}

	noRare 		= FALSE
	if (is.null(rarepoints)){
		rarepoints 	= seq(200,1000,len=5);
	} else if (all(is.na(rarepoints))){
		noRare		= TRUE
	}
	rarepoints=sort(unique(rarepoints))
	if (length(samples)!=1){
		stop("\"samples\" has to be a single number");
	}
	if (samples>1){
		samples 	= 1:samples;
	} else {
		samples 	= 1;
	}
	M = getMatrix(Mat,type="Feature",col.subset=subset,poolSamples=poolSamples)

	#check if subfraction really exist
	sfc = subFeatures[subFeatures%out%"all"]
	if(length(sfc > 0 )){
		for (i in 1:length(sfc)){
			if(sum(grep(sfc[i],dimnames(M)[[1]])) < 1){
				if (!silent){
					writeInfo(paste("Could not find subFeature",sfc[i]),opt,"warning");
				}else{
					warning(paste("Could not find subFeature",sfc[i]));
				}
			}
		}
	}

	matWarn(M,lowerThan="Norm");

	writeInfo(paste("Parameters: method=",paste(method,collapse=",")," ; samples=",paste(samples,collapse=", ")," ; rarepoints=",paste(rarepoints,collapse=",")),Opt,write=!silent)

	adds = FALSE;

	if (!is.null(MetaInfo) && !is.null(poolSamples)){
		colObj 	= getSampleColorName(poolSamples,MetaInfo)
		colsB 	= cols=colObj$ColLegend
		cls 	= names(colObj$ColLegend);  names(cls) = cls; cls= factor(cls)
		#attr(colObj$ColLegend,"pch")
		ltys 	= attr(colObj$ColLegend,"lty")
		lcls 	= names(cols);	llcls = length(lcls);
	} else if (!is.null(MetaInfo)){
		#browser()
		colObj = getSampleColorName(M,MetaInfo,givenCol=cls)

		if (!all(is.na(colObj$SampleColorNames))){
			cls = colObj$SampleColorNames;
			#smplNames = colObj$SampleNames
			if(is.null(smplcol)) {smplcol = colObj$SampleColor}
			colsB 	= cols=colObj$ColLegend
			lcls 	= names(cols);	llcls = length(lcls);
			if (!is.null(Opt$report$grayScale) && Opt$report$grayScale){ #line style differs only in b/w
				ltys = attr(cols,"lty")
			} else {
				ltys = array(1,llcls)
			}
		}
	}
	if (is.null(cls)){
		cls = factor(array(1,dim(M)[2]))
		names(cls) = dimnames(M)[[2]]
		colsB=cols=array("black")
		llcls=ltys=1; lcls = levels(cls)
#	} else {
#
#		if (any(is.na(cls))){
#			cls = factor(cls[!is.na(cls)])
#		}
#
#		if (!is.factor(cls)){cls = factor(cls)}
#		lcls = unique(cls);	llcls = length(lcls);
#
#		colObj = getSampleColorName(cls,MetaInfo,givenCol)
#		if(is.null(smplcol)){smplcol = colObj$SampleColor}
#		cols = makeClassCol(cls,smplcol,lcls=lcls);

#		cols = getGeneralColors(llcls,Opt);
#		colsB = getGeneralColors(llcls,Opt,alpha=0.5);
#		names(cols) = names(colsB) = lcls;
#		ltys = array(1,llcls)
	}


	#subselections = .extractComplexFeature(subFeatures,dimnames(M)[[1]])
	for (sf in 1:length(subFeatures)){
		if (is.null(poolSamples)){
			M = getMatrix(Mat,"Feature",col.subset=names(cls),row.subset=subFeatures[sf])
		} else {
			#M=M
		}

		if (removeSubRareSmpls){
			rarepoints2 =rarepoints
			subs = colSums(M)>=max(rarepoints)
			M = M[,subs]
		}else{
			rarepoints2 = suitable_rarefactions(rarepoints,M,Opt,extendRseries,silent=silent);
		}
		#rarefaction function
		if (is.null(givenRichness)){
			if (is.null(poolSamples)){
				rr=Rarefaction.curve(M,Opt,size=rarepoints2,subset=subset,
					repeats = samples, method = method,Mat = Mat,
					subFraTag = subFeatures[sf],Cl=Cl,recalcRare=recalcRare,
					tree=PHtree,noRare=noRare);
			} else {
				rr=Rarefaction.curve(M,Opt,size=rarepoints2,subset=subset,
						repeats = samples, method = method,Mat = M,
						subFraTag = subFeatures[sf],Cl=Cl,recalcRare=recalcRare,
						tree=PHtree,noRare=noRare);

			}
		} else {
			givenRichness = givenRichness[!is.na(givenRichness)]
			rr=rbind(givenRichness);
		}
		#main plot with sub sample dependent rarefaction curve
#			plotTy=plotType#single,box,vio
		#writeInfo(paste("Rarefaction on",subFeatures[sf]),opt)
		if (dim(rr)[[1]] == 1){
			maxIdx=1
		} else{
			nonNAs = apply(!is.na(rr),1,sum)
			maxIdx = max(which(max(nonNAs)==nonNAs))
		}
		if (subFeatures[sf]=="all"){
			#lastEst = rr[dim(rr)[1],];
			#ret[["Counts"]] = lastEst;
			ret[["Estimate"]] = rr;
			ret[["MaxEstimate"]] = rr[maxIdx,]
			ret[["Depth"]] = rarepoints2;
		} else {
			r2 = list();
			#lastEst = rr[dim(rr)[1],];
			#r2[["Counts"]] = lastEst;
			r2[["Estimate"]] = rr;
			r2[["MaxEstimate"]] = rr[maxIdx,]
			r2[["Depth"]] = rarepoints2;
			ret[["subfeatures"]][[subFeatures[sf]]] = r2;
			if (is.null(ret[["Depth"]])){
				ret[["Depth"]] = rarepoints2;
			}
		}
	}
	m_title = paste(whatIs,attr(rr,"methodType"))

	.plot.rarefaction(ret[["Estimate"]],cls,ret[["Depth"]],Opt,MetaInfo,
			cols=cols,plotTy=plotType,silent=silent,whatIs=whatIs,ltys=ltys,
			compactLeg = TRUE,poolSamples=poolSamples)

#		if (!RarePlot){
#			m_title = paste(whatIs,"diversity estimate")
#		}

	#plot density of highest samples
	if (!silent && maxEstPlot=="density"){
		.plot.rarefaction.density(ret,cls,ret[["Depth"]],Opt,MetaInfo,
				cols,colsB,silent,m_title=m_title)
	} else if (!silent && maxEstPlot=="box"){
		pix = c(400,400)
		writeInfo(paste("Boxplot of class",attr(rr,"methodType"),"at",ret[["Depth"]][maxIdx],"rarefied sequences"),Opt,"pic",second=pix)
		bpl = boxSingleWin(ret[["MaxEstimate"]],Opt,factor(cls),MetaInfo,leg.ty="below",
				OrderByMedian=T,main=m_title,smplCol=smplcol,block=blck[names(cls)])
		writeInfo("",Opt,"dev.off")
		aveSdClsMat(ret[["MaxEstimate"]],factor(cls),Opt,bpl$lcls,perc=F)
	}

	#writeInfo("Taxa Richness estimation",Opt,"subsectionfile");
	ret[["classes"]] = lcls

	if (DOtest && llcls > 1){
		M = ret[["MaxEstimate"]];

		#writeInfo("Testing mean and Variation",Opt,"subsectionfile")
		KT = .Univar.Core(M,class_=factor(cls),Opt=Opt,method="kruskal",
				block_=blck[names(cls)])
		KT$par$title="Test for equal means";
		LT = .Univar.Core(M,class_=factor(cls),Opt=Opt,method="levene")
		LT$par$title="Test for equal variance";
		if (!silent){
			Univar.Output.Singl(KT,Opt)
			Univar.Output.Singl(LT,Opt)
		}
		#writeInfo("",Opt,"mainfile")
	}

	#if (!silent){writeInfo("",Opt,"mainFile",write=subsection);}
	return(invisible(ret));
}
