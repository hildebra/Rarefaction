#include "Rwrapper.h"
//#include "IO.h"
#include "Rare.h"




#include <Rcpp.h>
using namespace Rcpp;



// get the options passed from R
//
// parameters
// mat : matrix object or path to file


// helper function, that just converts a div
// into an R list, as we have to do this serveral times
List createDivList(DivEsts * div){
    List divLst         = List::create(
		                      Named("samplename",   div->SampleName),
		                      Named("richness",     div->richness),
		                      Named("shannon",      div->shannon),
		                      Named("simpson",      div->simpson),
		                      Named("invsimpson",   div->invsimpson),
		                      Named("chao1",        div->chao1),
		                      Named("eve",          div->eve));
    return divLst;
}

IntegerMatrix matrix2Mat(std::vector<vector<unsigned int>> dfMat,
						std::vector<string> colnames, std::vector<string> rownames ){
	// create a mat from a vector vector uint
	IntegerMatrix NM 	= Rcpp::IntegerMatrix( dfMat[0].size(), dfMat.size());
	for (int i = 0; i < NM.ncol(); i++) {
		Rcpp::NumericVector  v = wrap(dfMat[i]);
		NM(_,i) = v;
	}

	Rcpp::List dimnms =  Rcpp::List::create(
						rownames, colnames);

	NM.attr("dimnames") = dimnms;

	return(NM);
}


// here we define what R passes on and returns to the C++ code
// [[Rcpp::export]]
List rcpp_rarefaction(Rcpp::String input, Rcpp::String output,
                    NumericMatrix rMatrix, StringVector inColNames,
					StringVector inRowNames,
                    int repeats, long depth, int NoOfMatrices,
                    bool verbose = false, bool returnObject = false, int margin=2) {

    // initialize variables
    std::vector< std::vector < double > > rmat;
    vector < string > incolnames;
	vector < string > inrownames;

    if(input == ""){
        // use R matrix as input
        int nc = rMatrix.ncol();
        rmat.resize(nc);
        for( int i=0; i<nc; i++){
            NumericMatrix::Column col = rMatrix(_,i) ;
            rmat[i].assign( col.begin() , col.end() ) ;
        }
        // also assign the colnames, so we convert the type
        incolnames =  Rcpp::as<vector < string > >(inColNames);
		inrownames =  Rcpp::as<vector < string > >(inRowNames);
    }

	// transpose matrix, yes or no
	bool transpose = false;
	if(margin == 1){ 	// apply function over rows, which is not the default, but columns
						// is what the function currently dioes by default
		transpose = true;
	}

    // create variables to be filled
    DivEsts * dd(NULL);
    dd 				= new DivEsts();
    vector<DivEsts*> * divvs;
    divvs 			=  new vector<DivEsts*>;

    // return vector for counts
    //std::vector<vector<unsigned int>> retCnts;
	std::vector<vector<vector<unsigned int>>> retCnts(NoOfMatrices); // initialize a vector of matrices with the number of repeats
	std::vector<string> retCntsSampleNames;
	std::vector<string> rowNames;

    // call the rarefaction main function
    int res = rarefyMain(input, output, "rare_inmat", repeats, depth, verbose,
                        returnObject, rmat, incolnames, inrownames , dd, divvs, retCnts, retCntsSampleNames,
						rowNames, NoOfMatrices, transpose);


	if(verbose == true){
		cout << "Done rarefying, will now produce R objects in Cpp\n";
		cout << "and pass them to R\n";
	}

    // convert output to R
    List divLst = createDivList(dd);

    // list of all divs of each sample
    std::list<Rcpp::List> majorLst;
	//DataFrame DFdivvs(divvs->size());
	std::vector<Rcpp::IntegerMatrix> RrarefyMatrices(NoOfMatrices); // vector to hold te matrices
    if(returnObject == true){
      // only create these objects if needed
      // all divvs
      int i;
      for(i=0; i<divvs->size(); i++){
        // create a Lst from div pointer
        List tmpDivLst = createDivList((*divvs)[i]);
        majorLst.push_back(tmpDivLst);
		//DFdivvs[i] = createDivVec((*divvs)[i]);
      }

      // matrices with all the counts
	  for(i=0; i < NoOfMatrices; i++){
		  IntegerMatrix RdfTmp 		= matrix2Mat(retCnts[i], retCntsSampleNames, rowNames);
		  RrarefyMatrices[i]	= RdfTmp;
	  }
    }
	if(verbose == true){
		cout << "All R objects were produced\n";

	}
    // create R object to return to R
    List returnList;
	if(returnObject == true){
		List retMatDF = wrap(RrarefyMatrices);

		returnList            = List::create(Named("div", divLst),
                                            Named("divvs", majorLst),
											Named("raremat",retMatDF));
      // add more objects here
    }else{
      returnList            = List::create(Named("div", divLst));
    }

    // clean up variables
    delete dd;

    return returnList;
}
