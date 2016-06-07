path <- '/Users/saary/testData/testdataSingle.csv'
#path <- '/Users/saary/testData/testdataSingle.csv'

# benchmark the ram usage of rarefaction using vagrind
cat("vegan")
args = commandArgs(trailingOnly=TRUE)
cat(args[1])

require("vegan")
rrarefy <-
    function(x, sample, margin=1)
{
    if (!identical(all.equal(x, round(x)), TRUE))
        stop("function is meaningful only for integers (counts)")
	if(margin != 1 && margin != 2){
		warning("Margin must be one or two. Will be set to 1 now (default)")
		margin <- 1
	}
	x <- as.matrix(x)
    if (ncol(x) == 1)
        x <- t(x)
    if (length(sample) > 1 && length(sample) != nrow(x))
        stop(gettextf(
             "length of 'sample' and number of rows of 'x' do not match"))
	if(margin == 1){
		sample <- rep(sample, length=nrow(x))
		colnames(x) <- colnames(x, do.NULL = FALSE)
		nm <- colnames(x)
	    ## warn if something cannot be rarefied
	    if (any(rowSums(x) < sample))
	        warning("Some row sums < 'sample' and are not rarefied")
	    for (i in 1:nrow(x)) {
	        if (sum(x[i,]) <= sample[i]) ## nothing to rarefy: take all
	            next
	        row <- sample(rep(nm, times=x[i,]), sample[i])
	        row <- table(row)
	        ind <- names(row)
	        x[i,] <- 0
	        x[i,ind] <- row
	    }
	}else if(margin == 2){
    cat(c('margin=2'))
		sample <- rep(sample, length=ncol(x))
		rownames(x) <- rownames(x, do.NULL = FALSE)
		nm <- rownames(x)
	    ## warn if something cannot be rarefied
	    if (any(colSums(x) < sample))
	        warning("Some colum sums < 'sample' and are not rarefied")
	    for (i in 1:ncol(x)) {
	        if (sum(x[,i]) <= sample[i]) ## nothing to rarefy: take all
	            next
	        col <- sample(rep(nm, times=x[,i]), sample[i])
	        col <- table(col)
	        ind <- names(col)
	        x[,i] <- 0
	        x[ind, i] <- col
	    }
	}
	#sample <- rep(sample, length=nrow(x))


    x
}

# load data in matrix
data     		<- read.table(file = path, sep = "\t", header = T, row.names = 1)
#samplesize 		<- min(colSums(data))

result.rarefy   <- rrarefy(x = data, sample = args[1], margin = 2)
str(result.rarefy)
