r.median <- function(n, x){
  median(unlist(sapply(x, function(x){x[[n]]}), F, F))
}

rare.status <- function(msg, verbose=TRUE){
  if(verbose == TRUE){
    cat(msg)
  }
}
# this is the rarefy function
rare <- function(input, repeats = 10, depth = 0, ReturnMatrix = 0, margin = 2, verbose = FALSE, threads = 1, tmpdir = NULL ){


	if(repeats < ReturnMatrix ){
		repeats <- ReturnMatrix
		warning(paste("Repeats can not be smaller than number of matrices to return. Repeats set to match ReturnMatrix. repeats = ReturnMatrix =", repeats, sep=" "))
	}

	#if(!all(depth > 0)){
	#	warning("Will now rarefy to 0.95 times the smallest column sum")
	#}

	# sort depths
	depth <- sort(as.numeric(depth))


	# convert dataframes
	if(class(input) == "data.frame"){
		input <- as.matrix(input)
	}

	#validate if input is a path or a matrix
	if(class(input) == "matrix"){
		rare.status("Matrix object supplied for analysis", verbose)
		# validate that the matrix is numeric
		if(!is.numeric(input)){
			stop("The supplied matrix object is not numeric. Please check your input matrix.")
		}

		# pass 1:x to Cpp as colnames
		removeCnames <- FALSE
		removeRnames <- FALSE
		if(is.null(colnames(input))){
			colnames(input) <- paste("col ", seq(1:ncol(input)), sep="")
			removeCnames <- TRUE
		}
		if(is.null(rownames(input))){
			rownames(input) <- paste("row ", seq(1:nrow(input)), sep="")
			removeRnames <- TRUE
		}

		# call the actual software

		result <- lapply(depth, function(d){
						res <- rcpp_rarefaction("", input, colnames(input),
								rownames(input),
			                    repeats, d,
								ReturnMatrix,
			                    verbose, threads,
								margin, "NULL", FALSE)

						# remove col and/or row names, as we've added them for
						# the Cpp software to work well
						if(removeRnames == TRUE && removeRnames == TRUE){
							res$raremat <- lapply(res$raremat, unname)
						}else{
							if(removeRnames == TRUE){
								res$raremat <- lapply(res$raremat, function(x) {rownames(x) <- NULL; return (x)})
							}
							if(removeCnames == TRUE){
								res$raremat <- lapply(res$raremat, function(x) {colnames(x) <- NULL; return (x)})
							}
						}
						return(res)
					})


	}else if(class(input) == "character"){
	  rare.status("A path to a matrix file was supplied", verbose)
    if(!is.null(tmpdir) & margin == 2){
      uselowmem <- TRUE;
      rare.status("Low memory mode will be used. Temporary files will be stored on storage medium", verbose)
    }else if(!is.null(tmpdir) & margin != 2){
      warning("Can not use low mem on margin = 1. Please consider transforming your input data, to use low mem mode.")
    }else{
      uselowmem   <- FALSE;
      tmpdir     <- "NULL";
    }

		# validate that the file exists
		if(!file.exists(input)){
			stop(paste("The file can not be found. Please verify that the file exists in the given location. The path given is:", input, sep = " "))
		}
		result <- lapply(depth, function(d){
							res <- rcpp_rarefaction( input,
	                        matrix(1,1,c(1)),
							c(NA),c(NA), # col and rownames
	                        repeats, d,
							ReturnMatrix,
	                        verbose, threads,
							margin, tmpdir, uselowmem)
							return(res)
						})
  }else{
    stop("Unknown input type. Path to a file (character) or a numeric matrix are accepted types.")
  }
  	result <- lapply(result, function(res){

		if(length(res$skipped) > 0){
			warning(paste(length(res$skipped), "samples where skipped because the depth was greater than the number of elements in the sample."))
		}

	  # calculate median for diversity measures
	  measures 					<- c('richness', 'shannon', 'simpson', 'invsimpson', 'chao1', 'eveness')
	  res$div.median 			<- lapply(measures, r.median, x=res$divvs)
	  names(res$div.median) 	<- paste("median.", measures, sep = "")

	  return(res)
	})

	if(length(depth) == 1){
		result <- result[[1]]
	}

	result$depths <- depth
	result$repeats <- repeats
	# set our class
	class(result) <- "rarefaction";
	return(result)

}
