# this is the rarefy function

rare <- function(input, output = ".", repeats=10, rareDepth = 1000,
				NoOfMatrices = 1,
				verbose=TRUE, returnObject=FALSE, margin = 2){

    # empty return object
    result <- list()

  #validate if input is a path or a matrix
  if(class(input) == "matrix"){
    rare.status("Matrix object supplied for analysis", verbose)

    # validate that the matrix is numeric
    if(!is.numeric(input)){
      warning("The supplied matrix object is not numeric. Please check your input matrix.")
    }

	# pass 1:x to Cpp as colnames
	if(is.null(colnames(input))){
		colnames(input) <- paste("col ", seq(1:ncol(input)), sep="")
	}
	if(is.null(rownames(input))){
		rownames(input) <- paste("row ", seq(1:nrow(input)), sep="")
	}

    result <- rcpp_rarefaction("", output,
                        input, colnames(input),
						rownames(input),
                        repeats, rareDepth,
						NoOfMatrices,
                        verbose, returnObject,
						margin)
  }else if(class(input) == "character"){
    rare.status("A path to a matrix file was supplied", verbose)

    # validate that the file exists
    if(!file.exists(input)){
      warning("The file can not be found.")
    }
    result <- rcpp_rarefaction(   input, output,
                        matrix(1,1,c(1)),
						c(NA),c(NA), # col and rownames
                        repeats, rareDepth,
						NoOfMatrices,
                        verbose, returnObject,
						margin)


  }else{
    warning("Unknown input type. Path to a file (character) or a numeric matrix are accepted types.")
  }

  return(result)

}

rare.status <- function(msg, verbose=TRUE){
  if(verbose == TRUE){
    cat(msg)
  }
}
