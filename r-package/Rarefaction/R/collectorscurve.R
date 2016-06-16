collectorscurveold <- function(x,  ...){
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


collectors.curve <- function(x, col = 0, times = 10, ...){
  if(class(x) == "rarefaction"){
    cum.sample.rare(x, col, times, ...)
  }else{
    cum.sample(x, col, times, ...)
  }
  return()
}



cum.sample.rare <- function(x, col = 1, times = 10,  ...){
  if(class(x) != "rarefaction"){
    stop("Not a rarefaction object")
  }
  
  depths <- x$depths
  
  if(length(depths) > 1){
    warning("might get cluttered")
  }
  
  a <- lapply(seq(1, length(depths), by=1), function(i, x){
    if(length(depths) == 1){
      b <- cum.sample(x$raremat[[1]], times = times, do.plot  = F)
    }else{
      b <- cum.sample(x[[i]]$raremat[[1]],  times = times, do.plot  = F)  
    }
    
    
    return(b)
  }, x = x)
  
  
  ymax <- max(unlist(a))
  ymin <- min(unlist(a))
  col <- rep_len(col, length(a))
  
  plot(1,type="n", , ylim=c(ymin,ymax), xlim=c(1, ncol(a[[1]])), ...) 
  mapply(function(ab, color){
    boxplot(ab,  add = TRUE, axes = FALSE, col = color)
  }, ab = a, color = col  )
  
  return()
}


cum.sample <- function(df, col= 0, times = 10, do.plot = TRUE, ...){ 
  # plot the diversity when picking random samples
  s <- seq(1, ncol(df), 1)
  m <- t(sapply(seq(1, times, 1), sample.diversity, df, s))
  if(do.plot ){
    boxplot(m, col = col , ...)  
  }
  rownames(m) <- seq(1, times, 1)
  colnames(m) <- s
  return(m)
}

sample.diversity <- function(i,df, s){
  d <- sapply(s,function(j){
    n <- sample(seq(1, ncol(df)), j)
    if(j >1){
      m <- rowSums(df[,n])
    }else{
      m <- n
    }
    length(m[m!=0]) 
  })
  return(d)
}




