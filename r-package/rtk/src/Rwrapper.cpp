#include "Rwrapper.h"
//#include "IO.h"
#include "RRare.h"




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
		                      Named("eveness",      div->eve));
    return divLst;
}

IntegerMatrix matrix2Mat(std::vector<rare_map>& dfMat,
						std::vector<string> colnames, std::vector<string> rownames, bool transpose=false ){
	// create a mat from a vector vector uint
	IntegerMatrix NM;
	Rcpp::List dimnms;
	if(transpose == false){
		NM 	= Rcpp::IntegerMatrix( rownames.size(), colnames.size());
		for (int i = 0; i < NM.nrow(); i++) {
      for(int j = 0; j < NM.ncol(); j++){
        auto fnd = dfMat[j].find(i);
        if(fnd != dfMat[j].end()){
          NM(i,j) = fnd->second;
        }else{
          NM(i,j) = 0;
        }
      }
		}
		dimnms =  Rcpp::List::create(rownames, colnames);
	}else{
    NM 	= Rcpp::IntegerMatrix( colnames.size(), rownames.size());
    for (int j = 0; j < NM.nrow(); j++) {
      for(int i = 0; i < NM.ncol(); i++){
        auto fnd = dfMat[j].find(i);
        if(fnd != dfMat[j].end()){
          NM(j,i) = fnd->second;
        }else{
          NM(j,i) = 0;
        }
      }
    }
		dimnms =  Rcpp::List::create( colnames, rownames);
	}

	// assign dimnames
	NM.attr("dimnames") = dimnms;

	return(NM);
}


// here we define what R passes on and returns to the C++ code
// [[Rcpp::export]]
List rcpp_rarefaction(Rcpp::String input,
						NumericMatrix rMatrix, StringVector inColNames,
						StringVector inRowNames,
						int repeats, long depth, int NoOfMatrices,
						bool verbose = false, unsigned int threads = 1,
						 int margin=2, Rcpp::String tmpDir = "", bool lowmem = false)
						{

	// check for user interrup
	Rcpp::checkUserInterrupt();

	// initialize variables
	std::vector< std::vector < double > > rmat;
	vector < string > incolnames;
	vector < string > inrownames;
  string outF = "";
  string mode = "rare_inmat";

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

  // switch to low mem if wanted
  if(lowmem == true){
    mode = "rare_lowMem";
    outF = tmpDir;
  }else{
    mode = "rare_inmat";
  }

	// transpose matrix, yes or no
	bool transpose = false;
	if(margin == 1){ 	// apply function over rows, which is not the default, but columns
						// is what the function currently dioes by default
		transpose = true;
	}

    // create variables to be filled
    vector<DivEsts*>  divvs(0,NULL);
	// return vector for counts
	std::vector<vector<rare_map>> retCnts(NoOfMatrices); // initialize a vector of matrices with the number of repeats
	std::vector<string> retCntsSampleNames;
	std::vector<string> rowNames;
	std::vector<string> skippedSamples;
  std::vector<double> ACE;
  vector<mat_fl> ICE;
  std::vector<double> chao2;

	// check for user interrup
	Rcpp::checkUserInterrupt();

	// call the rarefaction main function
	rarefyMain(input, outF, mode,  repeats, depth,  threads, verbose,
				 rmat, incolnames, inrownames ,
				 divvs, retCnts, retCntsSampleNames,
				 skippedSamples, ACE, ICE, chao2,
				rowNames, NoOfMatrices, transpose);

	// check for user interrup
	Rcpp::checkUserInterrupt();
	if(verbose == true){
		Rcout << "Done rarefying, will now produce R objects in Cpp\n";
		Rcout << "and pass them to R\n";
	}

    // convert output to R
    //List divLst = createDivList(dd);

    // list of all divs of each sample
    std::list<Rcpp::List> majorLst;
	if(verbose == true){
		Rcout << "Will now prepare diversity measures for R\n";
	}
	for(uint i = 0; i < divvs.size(); i++){
		// create a Lst from div pointer
		List tmpDivLst = createDivList(divvs[i]);
		majorLst.push_back(tmpDivLst);
		//delete divvs[i];
	}



	std::vector<Rcpp::IntegerMatrix> RrarefyMatrices(NoOfMatrices); // vector to hold te matrices

    if(NoOfMatrices > 0){
      // matrices with all the counts
	  if(verbose == true){
		  Rcout << "Will now prepare rarefied matrices for R\n";
	  }
	  for(int i=0; i < NoOfMatrices; i++){
		  if(retCnts[i].size() > 0){
			  IntegerMatrix RdfTmp 	= matrix2Mat(retCnts[i], retCntsSampleNames, rowNames, transpose);
			  RrarefyMatrices[i]		= RdfTmp;
		  }
	  }
    }
	if(verbose == true){
		Rcout << "All R objects were produced\n";
	}

    // create R object to return to R
    List returnList;

	if(NoOfMatrices > 0 ){
		List retMatDF;
		retMatDF 			= wrap(RrarefyMatrices);
		returnList			= List::create(	Named("divvs", majorLst),
											Named("raremat",retMatDF),
                      Named("ICE", wrap(ICE)),
                      Named("ACE", wrap(ACE)),
                      Named("chao2", wrap(chao2)),
											Named("skipped", wrap(skippedSamples)));
    }else{
		returnList			= List::create(Named("divvs", majorLst),
                      Named("ICE", wrap(ICE)),
                      Named("ACE", wrap(ACE)),
                      Named("chao2", wrap(chao2)),
										  Named("skipped", wrap(skippedSamples)));
  }


    return returnList;
}
