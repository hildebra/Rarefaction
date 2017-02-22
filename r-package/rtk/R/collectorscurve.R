#' x=primary input matrix
#' y=secondary input matrix for comparative plots
#' @param col: colors used for rarefaction plots (borders of boxplots, or both fill and border if col2 is left unspecified)
#' @param col2: colors used for rarefaction plots (fill of boxplots)
#' @param bin: bin size for collectors curve
#' @param add: add the plot to an existing plot?
#' @param doPot: should this function plot the collectors curve, or just return an object that can be plotted later with this function?
#' @param rareD: if given, rarefy input matrix to given depth prior to calculating the collectors curve
#' @param cls: vector describing the class of each input sample in matrix y. Will split the rarefaction curves for each class or accumulate successively within each class, if accumOrder is given
#' @param accumOrder: accumulate successively within each class, given by cls in the order given in this vector. All classes in cls must be represented in this vector.
collectors.curve = function(x, y=NULL, col = 1, times = 10, bin = 3, add=FALSE, ylim=NULL, xlim=NULL, doPlot=TRUE, rareD=NULL, cls=NULL, pch=20, col2=NULL,accumOrder=NULL, ...){
    
    # check input for the case of accumulation Order curve
    if (!is.null(accumOrder)){
        if (is.null(cls)){
            stop("accumOrder argument also requires cls argument with matching items")
        }
    }
    # set colors2 equal to colors1 as we seem to have them
    if (is.null(col2) && !is.null(col)){
        col2=col
    }
    # get the dimension names for later
    if (class(x) != "colCurve" && !is.null(cls) ){
        if (!is.null(names(cls))){
            cls = cls[dimnames(x)[[2]]]
        } else if ( length(cls) == dim(x)[2] ) {
            names(cls) <- dimnames(x)[[2]]
        }
    }
    
#    if(class(x) == 'rtk'){
 #       if(length(x$depths) == 1){
  #          a <- x$raremat
   #         x <- a
    #        str(x)
     #   }
    #}

    # in case of colCurveObject, we can reuse older stuff
    if (class(x) == "colCurve"){
        cat("collectors curve object provided\n")
        cls         <- attr(x,"cls")
        accumOrder  <- attr(x,"accumOrder")
        col         <- attr(x,"col")
        pch         <- attr(x,"pch")
        col2        <- attr(x,"col2")
        bin         <- attr(x,"bin")
        times       <- attr(x,"times");
        
    } else if (!is.null(rareD)){ #rarefaction required?
        
        cat(paste("Rarefying input matrix to",rareD))
        
        if (length(rareD)>1){
            # Question: two depths and y is null: What now?
            if (!is.null(y)){
                stop("2D collectors curve only works at a single rarefaction depth.\n")
            }
            # seems like we will have a complex return object
            x   <- rtk(x,1,rareD,TRUE)
        } else {
            x   <- rtk(x,1,rareD,TRUE)$raremat[[1]];
            if (!is.null(y)){
                y <- rtk(y,1,rareD,TRUE)$raremat[[1]];
            }
            if (!is.null(cls) && !is.null(dimnames(x))){
                cls <- cls[dimnames(x)[[2]]]
            }
        }
        cat("Done\n")
    }
    # in case of colCurve object to nothing
    if (class(x) == "colCurve"){
        a   <- x
        x   <- NULL
        gc()
    }else if (!is.null(y)){#double core mode
        yl  <- cum.sample.matched(x,y,bin,times,accumOrder)
        a   <- yl$d
        b   <- yl$d2
    } else if (class(x) == "rtk") {
        a = cum.sample.rare(x, col, times, bin, cls, accumOrder)
        # in case of two depth, we get a list of matrices
    }else {
        a = cum.sample(x, col, times, bin, cls,  accumOrder)
    }

    if (class(x) != "colCurve"){
        # Define all the attributes of the object colCurve
        attr(a,"cls")           <- cls
        attr(a,"accumOrder")    <- accumOrder
        attr(a,"col")           <- col
        attr(a,"pch")           <- pch
        attr(a,"col2")          <- col2
        attr(a,"bin")           <- bin
        attr(a,"times")         <- times
        class(a)                <- "colCurve"
    }
  
    if (doPlot){ #relies entirely on "a" object
        #unified graph interface
        if (is.null(ylim)){
            ylim <- c(max(1,min(unlist(a), T, F)),max(unlist(a), T, F))
        }
        
        if (!is.null(y)){
            #2d plot
            if (is.null(xlim)){xlim = range(b) }
            if (!add){
                plot(1, type = "n", ylim =ylim, xlim =xlim,  ...)
            }
            
            boxplot(a, add = TRUE, axes = FALSE, col = col, border=col, at = colMeans(b),pch=pch) #as.numeric(dimnames(a)[[2]])
            
            
        } else {
         
            if (!add){
                plot(1, type = "n", ylim =ylim, xlim = c(1, max(unlist( lapply(a,ncol) )) * 
                                        bin), ...)
                if (!is.null(cls) || !is.null(accumOrder)){
                    ucls=unique(cls);ucls = ucls[!is.na(ucls)]
                    if(length(col)>1 && !is.null(names(col))){
                        if (all(names(col)%in%ucls)){col=col[ucls]}
                        legend("topleft",legend=names(col),fill=col)
                    } 
					if ( length(col2)>1 && !is.null(names(col2))){
                        if (all(names(col2)%in%ucls)){col2=col2[ucls]}
                        legend("topleft",legend=names(col2),fill=col2)
                    }
                }
            }
           
            if (!is.null(accumOrder)){
                if (length(col)>1 && !is.null(names(col))){col = col [accumOrder]}
                if (length(col2)>1 && !is.null(names(col2))){col2 = col2 [accumOrder]}
                if (is.null(col)){col=col2}
                if (is.null(col2)){col2=col}
                prev <- 1
                xx   <- 1
                for (i in 1:length(accumOrder)){
                    sumSel  <- sum(cls%in%accumOrder[i])
                    steps2  <- getStepsAccum(sumSel,bin)
                    sbset   <- prev:(prev+(length(steps2)-1))
                   
                    boxplot(a[[xx]][,sbset], add = TRUE, axes = FALSE, col = col[i], border=col2[i], at = jitter(as.numeric(colnames(a[[xx]]))[sbset],bin*0.5),pch=pch[i], boxwex= bin)
                    prev    <- max(sbset)+1
                }                    
                return(invisible(a))
            } else if (!is.null(cls)){
                if ( length(a)!=1 && length(col)==1){col <- rep_len(col[1], length(a))}
                if ( is.null(col2)){col2=col}
                if ( length(a)!=1 && length(col2)==1){col2 <- rep_len(col2[1], length(a))}
                if ( length(a)!=1 && length(pch)==1){pch <- rep_len(pch[1], length(a))}
                if ( is.null(col)){col=col2}
            }
            for (i in 1:length(a)){
                #str(as.numeric(colnames(a[[i]])))
                boxplot(a[[i]], add = TRUE, axes = FALSE, col = col[i], border=col2[i], at = jitter(as.numeric(colnames(a[[i]])),bin*0.5),pch=pch[i],boxwex= bin)
            }                    
        }
    }
    return(invisible(a))
}


getStepsAccum <- function(n,bin){
    # function to return binned steps until n
    steps <- seq(1, n, bin)
    if (steps[length(steps)] != n){
        steps <- c(steps,n)
    }
    return(steps)
}

sample.diversity <- function (i, df, bin = 5, cls = NULL,accumOrder=NULL){
    #cat(paste0(i," "));
    
    # initialise empty variables
    n       <- ncol(df)
    d       <- c()
    clsret  <- c();
    
    if(is.null(accumOrder)){
        sn    <- seq(1, n, 1)
        steps <- getStepsAccum(n, bin)
        d     <- c(d,sapply(steps, function(j, sn) {
                            if( length(sn)==1 || length(sn) == j ){
                                nx=sn
                            }else{
                                nx <- sample(sn, j, replace = FALSE, prob = NULL)
                            }
                            sum(rowSums(df[, nx,drop=FALSE]) != 0)
                        }, sn = sn))
    }else{
        tot    <- array(FALSE,dim(df)[[1]])
        steps  <- c()
        totsum <- 0
        
        for (u in 1:length(accumOrder)){
            sn      <-which(cls%in%accumOrder[u])
            steps2  <- getStepsAccum(length(sn),bin)

            d       <- c(d,sapply(steps2, function(j, sn) {
                                if( length(sn)==1 || length(sn)==j ){
                                    nx=sn
                                }else {
                                    nx <- sample(sn, j, replace = FALSE, prob = NULL)
                                }
                                sum( tot | (rowSums(df[, nx,drop=FALSE]) != 0) )
                            }, sn = sn))
                            
            tot     <- tot|apply(df[,sn,drop=FALSE]>0,1,any)
            steps   <- c(steps,steps2+totsum);
            totsum  <- totsum+length(sn)

        }
    }
    names(d) <- steps
    #cat(d);cat("\n")
    #attr(d,"clsn")=clsret
    return(d)
}
cum.sample.rare <- function (x, col = 1, times = 10, bin, cls=NULL, accumOrder=NULL){
    if (class(x) != "rtk") {
        stop("Not a rarefaction object")
    }
    depths <- x$depths
    if (length(depths) > 1) {
        stop("Only one rarefaction depth is supported")
        #warning("might get cluttered, because more than one matrix was supplied (more than one depth)")
    }
    a <- lapply(seq(1, length(depths), by = 1), function(i, x) {
                if (length(depths) == 1) {
                    b <- cum.sample(x$raremat[[1]], times = times, bin = bin, cls = cls,
                            do.plot = F)
                }
                else {
                    b <- cum.sample(x[[i]]$raremat[[1]], times = times, 
                            bin = bin, cls=cls, do.plot = F)
                }
                return(b)
            }, x = x)
    return(a[[1]])
}

cum.sample.matched = function( dx, dy, bin = 5, times=10,accumOrder=NULL){
    #browser()
    if (!is.null(accumOrder)){stop("accumOrder is not supported for 2D collector's curves")}
    n <- ncol(dx)
    sn <- seq(1, n, 1)
    steps <- seq(1, n, bin)
    r1 = r2 = matrix(0,0,length(steps)); 
    for (x in seq(1,times,1)){
        d           <- array(0,length(steps))
        names(d)    <- steps; d2= d;
        cnt = 1;
        for (j in steps){
            n <- sample(sn, j, replace = FALSE, prob = NULL)
            if (j > 1) {
                m <- rowSums(dx[, n])
                m2 <- rowSums(dy[, n])
            } else {
                m2 <- m <- n
            }
            d[cnt] = length(which(m != 0))
            d2[cnt] = length(which(m2 != 0))
            cnt = cnt+1
        }
        r1=rbind(r1,d);r2=rbind(r2,d2);
    }
    dimnames(r1)[[2]] = dimnames(r2)[[2]] = steps  
    return(list(d=r1,d2=r2))
}

cum.sample <- function (dfx, col = 0, times = 10, bin, cls, accumOrder=NULL, do.plot = TRUE){   
    useCls=FALSE
    if (!is.null(cls)){
        ucls <- unique(cls)
        ucls <- ucls[!is.na(ucls)]
        
        if (!is.null(accumOrder)){
            if (!all(accumOrder%in%ucls)){
                stop("accumOrder arguments need to be present in cls vector")
            }
            ucls   <- "all"
            useCls <- FALSE#mem saver 
        } else {
            useCls <- TRUE
        }
        
    } else {
        nc   <- ncol(dfx)
        cls  <- array("all",nc)
        ucls <- "all"
    
    }
    m=list();# 
    for (u in 1:length(ucls)){
        
        if (useCls){
            df = dfx[,which(cls==ucls[u])]            
        } else {
            df  <- dfx
            dfx <- NULL
        }
        
        gc();
        m [[ucls[u] ]] =t(sapply(seq(1, times, 1), sample.diversity, df, bin, cls,accumOrder))
    }
    return(m)
    
}




