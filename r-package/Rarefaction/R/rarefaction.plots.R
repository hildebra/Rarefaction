
plot.rarefaction <- function(x, div = c("richness"),  groups = NA, col = NULL, lty = 1, pch = NA, fit = "arrhenius", legend = TRUE, legend.pos = "topleft", log.dim = "", ...){
  
  if(!div%in% c('richness', 'shannon', 'simpson', 'invsimpson', 'chao1', 'eveness')){
    stop(paste("Not a possible plotting option for div:", div))
  }

  
  if(!fit%in%c("arrhenius", "michaelis-menten", "logis", FALSE)){
    warning(paste("Not a valid fitting model:", div,  "No fitting will be done."))
  }
  
	if(length(x$depth) == 1){
	  singlePlot(x,div, groups, col, lty, pch, legend, legend.pos,  ...)
	}else if(length(x$depth) > 1){
    # by default make rainbow colors
    if(is.null(col)){col <-  rainbow(length(x[[1]]$divvs))}
    # call the correct plot function for multiple rarefaction curves
	  multiPlot(x, div, groups, col, lty, pch, fit, legend, legend.pos, log.dim, ...)
	}else{
    warning("No depths provided")
	}
  return()
}


singlePlot <- function(obj, div, groups, col , lty, pch, legend, legend.pos , ...){
  
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
  usr <- par( "usr" )
  # kuskal test
  if(ncol(df) > 1){
    k <- kruskal.test(df)
    text((usr[2]-usr[1])/2, usr[4], adj = c(0.5,1.5), labels = paste("Kruskal-Wallis p-value = ", round(k$p.value,5)))
  }
  
  if(legend){
    legend(legend.pos, inset=.02,names(df) , fill=col, horiz=FALSE, cex=0.8)
  }
  
  
  return()
}




multiPlot <- function(obj, div, groups, col , lty, pch, fit,  legend, legend.pos, log.dim, ...){
	# takethe obj when there are more than one repeat
  depths <- obj$depths
  ydata <-  matrix(sapply(seq(1, length(depths), by=1), getDivvs, obj=obj, divName=div), ncol=length(depths), byrow = FALSE)
  colnames(ydata) <- paste("depth",depths)
  rownames(ydata) <- sapply(obj[[1]]$divvs, function(x){return(x$samplename)})
  
  # merge samples if grouping should be done
  if(!is.na(groups) && length(groups) == nrow(ydata)){
    ydata <- as.data.frame(ydata)
    
    ydata <- lapply((split(ydata, groups)), function(g){
        return(sapply(g, median, na.rm=TRUE))      
    })
    ydata <- as.data.frame(ydata)
    
  }else{
  ydata <- as.data.frame(t(ydata))
  }
  

  ymax <- max(unlist(ydata)[!is.na(unlist(ydata))])
  ymin <- min(unlist(ydata)[!is.na(unlist(ydata))])
  
  # set the pch and col to same length as ncol(ydata)
  col <- rep_len(col, ncol(ydata))
  pch <- rep_len(pch, ncol(ydata))  
  
  
  
  # initial plot
  plot(1,type="n",log = log.dim, ylim=c(ymin,ymax), xlim=c(min(depths),max(depths)), ...) # ,xlab=xlabel,ylab=ylabel,ylim=ylime,xlim=c(min(rarepoints),max(rarepoints)))
  usr <- par( "usr" )
  # plot the lines, one at a time
  legendcolors <- mapply(function(y, col, lty, pch, depths ){
    
    if(fit%in%c("arrhenius", "michaelis-menten", "logis")){
      
      # find the fitting model
      FITM <- function(fit, depths, y) {
        df <-  data.frame(x=depths, y = y)
        switch(fit,
          "arrhenius" = nls(y ~ fit.arrhenius(x,a,b), data = df, start = list(a = 1, b = 100)),
          "michaelis-menten"  = nls(y ~ SSmicmen(x, Vm, k), data = df),
          "logis" = nls(y ~ SSlogis(x, Asym, xmid, scal), data = df)
        )
      }
      
      # plot a smooth line with n = 2000 datapoints, experimanetally set number, might be valid as a parameter?
      xD <- seq(from = min(depths), to = max(depths), length.out = 2000)
      lines(xD, predict(FITM(fit, depths, y), list(x = xD)), col = col)
    }else{
      # no fitting
      lines(depths, y, lwd=1, col=col, lty = lty)  
    }

    if(!is.na(pch)){
      points(depths, y, col=col, lty = lty, pch = pch)  
    }
    
    return(col)
  }, y=ydata,col=col, lty, pch, MoreArgs=list(depths=depths))
  
  
  if(ncol(ydata) > 1){
    k <- kruskal.test(ydata[apply(ydata, 1, function(x){all(!is.na(x))}),])
    legend("top",  paste("Kruskal-Wallis p-value = ", round(k$p.value,5)), bty ="n", pch=NA)
  }
  
  if(legend){
    legend(legend.pos, inset=.02,names(ydata) , fill=legendcolors, horiz=FALSE, cex=0.8)
  }
  
  return(NULL)
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


fit.arrhenius <- function(x, a, b){
  v <- b * x^a
  return(v)
}



