collectors.curve <- function(x, col = 0, times = 10, bin = 5, ...){
  if(class(x) == "rtk"){
    cum.sample.rare(x, col, times, bin, ...)
  }else{
    cum.sample(x, col, times, bin, ...)
  }
  return()
}



cum.sample.rare <- function(x, col = 1, times = 10,  bin , ...){
  if(class(x) != "rtk"){
    stop("Not a rarefaction object")
  }
  depths <- x$depths

  if(length(depths) > 1){
    warning("might get cluttered, because more than one matrix was supplied (more than one depth)")
  }

  a <- lapply(seq(1, length(depths), by=1), function(i, x){
    if(length(depths) == 1){
      b <- cum.sample(x$raremat[[1]], times = times, bin = bin, do.plot  = F)
    }else{
      b <- cum.sample(x[[i]]$raremat[[1]],  times = times, bin= bin, do.plot  = F)
    }
    return(b)
  }, x = x)


  ymax <- max(unlist(a), T, F)
  ymin <- min(unlist(a), T, F)
  col <- rep_len(col, length(a))

  plot(1,type="n", , ylim=c(ymin,ymax), xlim=c(1, ncol(a[[1]])*bin), ...)
  mapply(function(ab, color){
    boxplot(ab,  add = TRUE, axes = FALSE, col = color, at = as.numeric(colnames(ab)))
  }, ab = a, color = col  )
  return()
}


cum.sample <- function(df, col= 0, times = 10, bin, do.plot = TRUE, ...){
  # plot the diversity when picking random samples
  #s <- seq(1, ncol(df), 1)
  m <- t(sapply(seq(1, times, 1), sample.diversity, df, bin))
  if(do.plot ){
    boxplot(m, col = col , ...)
  }
  #rownames(m) <- seq(1, times, 1)
  #colnames(m) <- s
  return(m)
}

sample.diversity <- function(i,df, bin = 5){
  n 		<- ncol(df)
  sn 		<- seq(1, n, 1)
  steps 	<- seq(1, n, bin)
  d <- sapply(steps,function(j, sn){

    n <- sample(sn, j, replace = FALSE, prob = NULL)
	if(j >1){
      m <- rowSums(df[,n])
    }else{
      m <- n
    }
    length(which(m != 0))
  }, sn = sn)
  names(d) <- steps
  return(d)
}
