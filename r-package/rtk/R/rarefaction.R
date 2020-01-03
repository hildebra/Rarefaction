r.median <- function(n, x){
  x1 <- unlist(sapply(x, function(x){x[[n]]}), F, F)
  x2 <- x1[!is.na(x1)]
  return(median(x2))
}

rare.status <- function(msg, verbose=TRUE){
  if(verbose == TRUE){
    cat(msg)
  }
}
# this is the rarefy function
rtk <- function(input, repeats = 10, depth = 1000, ReturnMatrix = 0, margin = 2, verbose = FALSE, threads = 1, tmpdir = NULL ){


    # pass 1:x to Cpp as colnames
    removeCnames <- FALSE
    removeRnames <- FALSE

  if(repeats < ReturnMatrix ){
    repeats <- ReturnMatrix
    warning(paste("Repeats can not be smaller than number of matrices to return. Repeats set to match ReturnMatrix. repeats = ReturnMatrix =", repeats, sep=" "))
  }


  # sort depths
  depth <- sort(as.numeric(depth))

  # convert dataframes
  if(is.data.frame(input)){
    input <- as.matrix(input)
  }

  #validate if input is a path or a matrix
  if(is.matrix(input)){
    rare.status("Matrix object supplied for analysis", verbose)
    # validate that the matrix is numeric
    if(!is.numeric(input)){
      stop("The supplied matrix object is not numeric. Please check your input matrix.")
    }
    if(any(is.na(input))){
        stop("The input data contains NA values. Please sanitize your input first.")
    }

    if(is.null(colnames(input))){
      colnames(input) <- paste("col ", seq(1:ncol(input)), sep="")
      removeCnames <- TRUE
    }
    if(is.null(rownames(input))){
      rownames(input) <- paste("row ", seq(1:nrow(input)), sep="")
      removeRnames <- TRUE
    }

    # call the actual software

    result <-  rcpp_rarefaction("", input, colnames(input),
                                    rownames(input), repeats, depth,
                                    ReturnMatrix, verbose, threads,
                                    margin, "NULL", FALSE)
    gc()

    
        
          


  }else if(is.character(input)){
    rare.status("A path to a matrix file was supplied", verbose)
    if(!is.null(tmpdir) & margin == 2){
      uselowmem <- TRUE;
      rare.status("Low memory mode will be used. Temporary files will be stored on storage medium", verbose)
    }else if(!is.null(tmpdir) & margin != 2){
      warning("Can not use low mem on margin = 1. Please consider transforming your input data, to use low memory (swap) mode.")
    }else{
      uselowmem   <- FALSE;
      tmpdir     <- "NULL";
    }

    # validate that the file exists
    if(!file.exists(input)){
      stop(paste("The file can not be found. Please verify that the file exists in the given location. The path given is:", input, sep = " "))
    }
    result <- rcpp_rarefaction( input,
                               matrix(1,1,c(1)),
                               c(NA),c(NA), # col and rownames
                               repeats, depth,
                               ReturnMatrix,
                               verbose, threads,
                               margin, tmpdir, uselowmem)

            
  }else{
    stop("Unknown input type. Path to a file (character) or a numeric matrix are accepted types.")
  }
    result <- lapply(result, function(res){
        # remove names, if there werent any
        if(removeRnames == TRUE && removeCnames == TRUE){
          res$raremat <- lapply(res$raremat, unname)
        }else{
          if(removeRnames == TRUE){
            res$raremat <- lapply(res$raremat, function(x) {rownames(x) <- NULL; return (x)})
          }
          if(removeCnames == TRUE){
            res$raremat <- lapply(res$raremat, function(x) {colnames(x) <- NULL; return (x)})
          }
        }
            
        

        if(length(res$skipped) > 0){
          warning(paste(length(res$skipped), "samples where skipped because the depth was greater than the number of elements in the sample."))
        }

        # calculate median for diversity measures
        measures               <- c('richness', 'shannon', 'simpson', 'invsimpson', 'chao1', 'eveness')
        res$div.median         <- lapply(measures, r.median, x=res$divvs)
        names(res$div.median)  <- paste("median.", measures, sep = "")

    return(res)
  })

  if(length(depth) == 1){
    result <- result[[1]]
  }else{
    names(result) <- depth
  }

  result$depths    <- depth
  result$repeats   <- repeats
  class(result)    <- "rtk";
  gc()
  return(result)
}






get.diversity <- function(obj, div = 'richness', multi = FALSE){
    if(multi == FALSE){
      if(!is(obj,'rtk')){
        stop("this function requires an object of type 'rtk'")
      }
    }
    if(length(obj$depths) == 1 || multi == T){
      if(multi == T){
        obj <- obj[[1]]
      }
      y <- sapply(obj$divvs, function(x){
          r <- x[[div]]
          return(r)
        })
    }else{
      y <- lapply(obj$depths, function(i){
          ia <- as.character(i)
          x <- get.diversity(obj[ia], div, multi = TRUE)
          return(x)
        })
    }
    return(y)
}

get.median.diversity <- function(obj, div = 'richness'){
   if(!is(obj,'rtk')){
      stop("this function requires an object of type 'rtk'")
    }
  res <- get.diversity(obj, div)
  if(length(obj$depths) == 1){
    ret <- apply(res, 2, median)
  }else{
    ret <- lapply(res, function(x){
        return(apply(x, 2, median))
      })
  }
  return(ret)
}
get.mean.diversity <- function(obj, div = 'richness'){
  if(!is(obj,'rtk')){
    stop("this function requires an object of type 'rtk'")
  }
  res <- get.diversity(obj, div)
  if(length(obj$depths) == 1){
    ret <- apply(res, 2, mean)
  }else{
    ret <- lapply(res, function(x){
        return(apply(x, 2, mean))
      })
  }
  return(ret)
}
