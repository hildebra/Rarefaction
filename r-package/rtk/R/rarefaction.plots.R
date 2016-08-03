
plot.rtk <- function(x, div = c("richness"),  groups = NA, col = NULL, lty = 1, pch = NA, fit = "arrhenius", legend = TRUE, legend.pos = "topleft", log.dim = "", boxplot = FALSE, ...){

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
    if(is.null(col)){
        if(!all(is.na(groups))){
          col <-  rainbow(length(unique(groups)))
        }else{
          col <-  rainbow(length(x[[1]]$divvs))
        }
    }
    # call the correct plot function for multiple rarefaction curves
    rarefaction.curve(x, div, groups, col, lty, pch, fit, legend, legend.pos, log.dim,boxplot,  ...)
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




rarefaction.curve <- function(obj, div, groups, col , lty, pch, fit,  legend, legend.pos, log.dim,boxplot, ...){
	# takethe obj when there are more than one depth
  depths <- obj$depths
  ydata <-  matrix(sapply(seq(1, length(depths), by=1), getDivvs.median, obj=obj, divName=div), ncol=length(depths), byrow = FALSE)
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


  ymax <- max(unlist(ydata, F,F)[!is.na(unlist(ydata,F,F))])
  ymin <- min(unlist(ydata, F,F)[!is.na(unlist(ydata,F,F))])

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
  }, y=ydata,col=col, lty, pch, MoreArgs=list(depths=sort(depths)))




  if(boxplot == TRUE){
    ydataB <-  lapply(seq(1, length(depths), by=1), getDivvs.raw, obj=obj, divName=div, rep=obj$repeats)
    # grouping here:
    groupsB <- groups
    if(all(is.na(groupsB))){
      # make groups for everything, so each group is size 1.
      groupsB <- 1:ncol(ydata)
    }
      a <- sapply(ydataB, function(y){

        y <- as.data.frame(t(y))
        y <- lapply(split(y, groupsB), function(g){
          return(apply(g,2, median, na.rm=TRUE))
        })

        return(y)
      })

      # plot the groupsB
    if(!is.null(dim(a))){
      a2 <- apply(a, 1, function(x){
        x <- as.data.frame(x)
        names(x) <- depths
        return(x)
      })
    }else{
      a2 <- list(a=as.data.frame(a))
    }
      mapply(function(x, color){
        boxplot(x, add= TRUE, fill = color, border = color, at = depths, axes=FALSE )
      }, x = a2, color = col)

  }



  if(ncol(ydata) > 1){
    k <- kruskal.test(ydata[apply(ydata, 1, function(x){all(!is.na(x))}),])
    legend("top",  paste("Kruskal-Wallis p-value = ", round(k$p.value,5)), bty ="n", pch=NA)
  }

  if(legend){
    legend(legend.pos, inset=.02,names(ydata) , fill=legendcolors, horiz=FALSE, cex=0.8)
  }

  return(NULL)
}






rarefaction.curve.boxplot <- function(x,  ...){
  if(class(x) != "rarefaction"){
    stop("Not a rarefaction object")
  }

  depths <- x$depths

  a <- sapply(seq(1, length(depths), by=1), function(i, ...){
    b <- sapply(x[[i]]$raremat, function(j){
      s <- rowSums(j)
      length(s[s != 0])
    })
    return(b)
  }, x = x)


  df <- as.data.frame(a)


  if(length(df) == 1){
    rownames(df) <- depths
    colnames(df) <- "V1"
    plot(depths, a , ...)
  }else if(length(df) > 1){
    colnames(df) <- depths
    boxplot(df, at = depths, ...)
  }

  return(df)
}












getDivMedians <- function(i, obj, divName){
  return(obj[[i]]$div.median[[paste("median.",divName, sep="")]])
}

getDivvs.median <- function(i, obj, divName){
  y <- sapply(obj[[i]]$divvs, function(x){
    median(x[[divName]])
  })
  return(y)
}
getDivvs.raw <- function(i, obj, divName, reps){
  y <- sapply(obj[[i]]$divvs, function(x){
    r <- x[[divName]]

    if(length(r) != reps){
      r <- rep_len(0, reps)
    }
    return(r)

  })
  return(y)
}


fit.arrhenius <- function(x, a, b){
  v <- b * x^a
  return(v)
}
