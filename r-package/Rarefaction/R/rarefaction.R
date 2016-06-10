r.median <- function(n, x){
  median(unlist(sapply(x, function(x){x[[n]]}), F, F))
}

rare.status <- function(msg, verbose=TRUE){
  if(verbose == TRUE){
    cat(msg)
  }
}
# this is the rarefy function
rare <- function(input, repeats=10, depth = 1000,
				NoOfMatrices = 0,
				margin = 2, verbose=TRUE, threads=1 ){

    # empty return object
    #result <- list()

	#output <- "." # no ouput

	if(repeats < NoOfMatrices){
		repeats <- NoOfMatrices
		warning(paste("Repeats can not be smaller than number of matrices to return. Repeats set to match NoOfMatrices. repeats = NoOfMatrices =", repeats, sep=" "))
	}

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
		result <- rcpp_rarefaction("", input, colnames(input),
							rownames(input),
		                    repeats, depth,
							NoOfMatrices,
		                    verbose, threads,
							margin)

		# remove col and/or row names, as we've added them for
		# the Cpp software to work well
		if(removeRnames == TRUE && removeRnames == TRUE){
			result$raremat <- lapply(result$raremat, unname)
		}else{
			if(removeRnames == TRUE){
				result$raremat <- lapply(result$raremat, function(x) {rownames(x) <- NULL; return (x)})
			}
			if(removeCnames == TRUE){
				result$raremat <- lapply(result$raremat, function(x) {colnames(x) <- NULL; return (x)})
			}
		}

	}else if(class(input) == "character"){
	  rare.status("A path to a matrix file was supplied", verbose)

		# validate that the file exists
		if(!file.exists(input)){
			stop(paste("The file can not be found. Please verify that the file exists in the given location. The path given is:", input, sep = " "))
		}
	    result <- rcpp_rarefaction( input,
	                        matrix(1,1,c(1)),
							c(NA),c(NA), # col and rownames
	                        repeats, depth,
							NoOfMatrices,
	                        verbose, threads,
							margin)
  }else{
    stop("Unknown input type. Path to a file (character) or a numeric matrix are accepted types.")
  }
  if(length(result$skipped) > 0){
	  warning(paste(length(result$skipped), "samples where skipped because the depth was greater than the number of elements in the sample."))
  }

  # calculate median for diversity measures
  measures 					<- c('richness', 'shannon', 'simpson', 'invsimpson', 'chao1', 'eve')
  result$div.median 		<- lapply(measures, r.median, x=result$divvs)
  names(result$div.median) 	<- paste("median.", measures, sep = "")


	# set our class
	class(result) <- "rarefaction";
	return(result)

}
