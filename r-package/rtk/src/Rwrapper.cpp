#include "Rwrapper.h"
//#include "IO.h"
#include "RRare.h"




#include <Rcpp.h>
using namespace Rcpp;



options::options(string in, string tmpDir, int r, std::vector<double> d, int NoOfMatrices, bool v, unsigned int t) :input(""), output(""), mode(""),
    referenceDir(""), referenceFile(""),
    depth(), repeats(10), write(0), threads(1), writeSwap(true), verbose(false),
    modDB(""), modRedund(5), modEnzCompl(0.5f), modModCompl(0.5f), modWrXtraInfo(false), 
    modCollapse(false), calcCoverage(false),
    xtra("")
    {
    depth  = d;
    input = in;
    output = tmpDir;
    repeats = r;
    write = NoOfMatrices;
    verbose = v;
    threads = t;
}

// get the options passed from R
//
// parameters
// mat : matrix object or path to file




// helper function, that just converts a div
// into an R list, as we have to do this serveral times
List createDivList(DivEsts * div, int di){
    List divLst         = List::create(
		                      Named("samplename",   div->SampleName),
		                      Named("richness",     div->richness[di]),
		                      Named("shannon",      div->shannon[di]),
		                      Named("simpson",      div->simpson[di]),
		                      Named("invsimpson",   div->invsimpson[di]),
		                      Named("chao1",        div->chao1[di]),
		                      Named("eveness",      div->eve[di]));
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



    
List returnRList(options* opts, vector<DivEsts*>& divvs, vector<vector<mat_fl>>& ACE, vector<vector<mat_fl>>& ICE, vector<vector<mat_fl>>& chao2, std::vector<string> skippedSamples, std::vector<string>  cntsNam, std::vector<string> rowNames, bool transpose, vector< vector< vector< rare_map > >> MaRare, unsigned int di){
    std::list<Rcpp::List> majorLst;
    
    for(uint i = 0; i < divvs.size(); i++){
	    // create a Lst from div pointer
	    List tmpDivLst = createDivList(divvs[i], di);
	    majorLst.push_back(tmpDivLst);
    }

	std::vector<Rcpp::IntegerMatrix> RrarefyMatrices(opts->write); // vector to hold te matrices

    if(opts->write > 0){
      // matrices with all the counts
	  if(verbose == true){
		  Rcout << "Will now prepare rarefied matrices for R\n";
	  }
	  for(int i=0; i < opts->write; i++){
		  if(MaRare[di][i].size() > 0){
			  IntegerMatrix RdfTmp 	= matrix2Mat(MaRare[di][i], cntsNam, rowNames, transpose);
			  RrarefyMatrices[i]		= RdfTmp;
		  }
	  }
    }

    // create R object to return to R
    List returnList;
	if(opts->write > 0 ){
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

// here we define what R passes on and returns to the C++ code
// [[Rcpp::export]]
List rcpp_rarefaction(Rcpp::String input,
						NumericMatrix rMatrix, StringVector inColNames,
						StringVector inRowNames,
						int repeats, Rcpp::NumericVector depth, int NoOfMatrices,
						bool verbose = false, unsigned int threads = 1,
						 int margin=2, Rcpp::String tmpDir = "", bool lowmem = false)
						{

	// check for user interrup
	Rcpp::checkUserInterrupt();
    // make options:
    options* opts = new options(input, tmpDir, repeats, Rcpp::as<std::vector<double>>(depth), NoOfMatrices, verbose, threads);

    
	// initialize variables
	std::vector< std::vector < double > > rmat;
	vector < string > incolnames;
	vector < string > inrownames;
    string outF = "";
    string mode = "memory";

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
    mode = "swap";
    outF = tmpDir;
    opts->writeSwap = true;
  }else{
    mode = "memory";
    opts->writeSwap = false;
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
	 vector< vector< vector< rare_map > >> MaRare(opts->depth.size(), vector< vector< rare_map> > (opts->write)); // initialize a vector of matrices with the number of repeats
    std::vector<vector<string>> cntsNames(opts->depth.size(), vector<string>());
	std::vector<string> rowNames;
	std::vector<string> skippedSamples;

	// store chao etc.
    vector<vector<mat_fl>> chao2(opts->depth.size());
    vector<vector<mat_fl>> ICE(opts->depth.size());
    vector<vector<mat_fl>> ACE(opts->depth.size());

	// check for user interrup
	Rcpp::checkUserInterrupt();
    if(verbose == true){
		Rcout << "\nPass data to C++ for rarefaction.\n";
		Rcout << "This might take long, depending on your input\n";
		Rcout << "please wait ...\n";
	}
	// call the rarefaction main function
	rarefyMain( opts, mode,
				 rmat, incolnames, inrownames ,
				 divvs, MaRare, cntsNames,
				 skippedSamples, ACE, ICE, chao2,
				rowNames, transpose);
						
	// check for user interrup
	Rcpp::checkUserInterrupt();
	if(verbose == true){
		Rcout << "\nDone rarefying, will now produce R objects\n";
		Rcout << "and pass them back to R\n";
	}

    // convert output to R
    if(verbose == true){
		Rcout << "Will now prepare diversity measures for R\n";
	}
	List returnList;
	for(unsigned int di = 0; di < opts->depth.size(); di++){
	    returnList.push_back(returnRList(opts, divvs, ACE, ICE, chao2, skippedSamples, cntsNames[di], rowNames, transpose, MaRare, di));
	}
    if(verbose == true){
		Rcout << "All R objects were produced\n";
	}

    return returnList;
}


options::options(int argc, char** argv) :input(""), output(""), mode(""),
    referenceDir(""), referenceFile(""),
    depth(), repeats(10), write(0), threads(1), writeSwap(true), verbose(false),
    modDB(""), modRedund(5), modEnzCompl(0.5f), modModCompl(0.5f), modWrXtraInfo(false), 
    modCollapse(false), calcCoverage(false),
    xtra("") {}
