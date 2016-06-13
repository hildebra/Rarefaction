
plot.rarefaction <- function(obj, div = c("richness"),  groups = NA, col = NULL, lty = 1, pch = NA, ...){
    
	if(length(obj$depth) == 1){
	  singlePlot(obj,div, groups, col, lty, pch, ...)
	}else if(length(obj$depth) > 1){
    # by default make rainbow colors
    if(is.null(col)){col <-  rainbow(length(obj[[1]]$divvs))}
    # call the correct plot function for multiple rarefaction curves
	  multiPlot(obj, div, groups, col, lty, pch, ...)
	}else{
    warning("No depths provided")
	}
}


singlePlot <- function(obj, div, groups, col , lty, pch, ...){
  
  # remove empty samples
  emptys <- sapply(obj$divvs, function(x){
    !x$samplename%in%obj$skipped
  })
  divvs <- obj$divvs[emptys]
  
  # get the data as df
  df <- sapply(divvs, function(x){
    return(unlist(x[[div]]))
  })
  colnames(df) <- sapply(divvs, function(x){
    return(unlist(x$samplename))
  })

  # grouping, is a little complex, but working fine I think
  if(!all(is.na(groups)) && (length(groups) == ncol(df) || length(groups) == length(emptys))){
    # remove unused groups
    if(length(groups) == length(emptys)){
      groups <- groups[emptys]
    }
    # make a dataframe for splitting
    df <- as.data.frame(t(df))
    # split and median
    df <- sapply(split(df, groups), function(gr){
      apply(gr, 1, median, na.rm=TRUE)
    })
  }else{
    df <- as.data.frame(df)
  }
  
  # set the pch and col to same length as ncol(ydata)
  if(is.null(col)){
    col <-  rainbow(length(df))
  }else{
    col <- rep_len(col, length(df))
  }
  
  pch <- rep_len(pch, length(df))  
    
  # boxplot the data
  boxplot(df, col=col, pch=pch, ...)

}




multiPlot <- function(obj, div, groups, col , lty, pch, ...){
	# takethe obj when there are more than one repeat
  depths <- obj$depths
  ydata <-  matrix(sapply(seq(1, length(depths), by=1), getDivvs, obj=obj, divName=div), ncol=length(depths), byrow = FALSE)
  colnames(ydata) <- paste("depth",depths)
  rownames(ydata) <- sapply(obj[[1]]$divvs, function(x){return(x$samplename)})
  
  # merge samples if grouping should be done
  if(!is.na(groups) && length(groups) == nrow(ydata)){
    ydata <- as.data.frame(cbind(ydata, groups))
    
    ydata <- lapply((split(ydata, ydata$groups)), function(g){
      sapply(g, median, na.rm=TRUE)
    })
    ydata <- as.data.frame(ydata)
    # remove groups again from DF object, so we have pure data
    ydata <- ydata[-which(rownames(ydata)=="groups"),]
    colnames(ydata) <- unique(groups)
  }else{
    ydata <- as.data.frame(t(ydata))
  }
  
  # get the y and x ranges
  ymax <- max(apply(ydata,1, function(x){max(x[!is.na(x)])}))
  ymin <- min(apply(ydata,1, function(x){min(x[!is.na(x)])}))
 
  # set the pch and col to same length as ncol(ydata)
  col <- rep_len(col, ncol(ydata))
  if(all(is.na(pch))){
    pch <- rep_len(pch, ncol(ydata))  
  }
  
  
  # initial plot
  plot(1,type="n",log="x", ylim=c(ymin,ymax), xlim=c(min(depths),max(depths)), ...) # ,xlab=xlabel,ylab=ylabel,ylim=ylime,xlim=c(min(rarepoints),max(rarepoints)))
  
  # plot the lines, one at a time
  legendcolors <- mapply(function(y, col, lty, pch, depths ){
  str(y)
    lines(depths, y, lwd=1, col=col, lty = lty)
    if(!is.na(pch)){
      points(depths, y, col=col, lty = lty, pch = pch)  
    }
    
    return(col)
  }, y=ydata,col=col, lty, pch, MoreArgs=list(depths=depths))

}



getDivMedians <- function(i, obj, divName){
  return(obj[[i]]$div.median[[paste("median.",divName, sep="")]])
}

getDivvs <- function(i, obj, divName){
  y <- sapply(obj[[i]]$divvs, function(x){
    median(x[[divName]])
  })
  return(y)
}
