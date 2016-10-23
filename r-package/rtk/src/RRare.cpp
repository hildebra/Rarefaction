
#include "RRare.h"


const char* rar_ver="0.63		 alpha R";

struct cDR{
	DivEsts* div;
	std::vector<vector<vector<uint>>> retCnts;
	string retCntsSampleName;
	vector<rare_map> RareSample;
	string skippedSample;
};

struct rareStruct{
	DivEsts* div;
	string cntsName;
	vector< rare_map> cnts;
	string skippedNames;
};

rareStruct* calcDivRar(int i, Matrix* Mo, DivEsts* div,  long rareDep,
	vector<vector<uint>>* abundInRow, vector<vector<uint>>* occuencesInRow,
	string outF,
	int repeats, int writeFiles,
	int NoOfMatrices){

	smplVec* cur 	= Mo->getSampleVec(i);
	string curS 	= Mo->getSampleName(i);
	div->SampleName = curS;

	// vector holding the rarefaction results for this sample
	// repeat times
	std::vector<rare_map> cnts;
	string cntsName;
	string skippedNames;
	cur->rarefy(rareDep, outF, repeats, div, cnts, cntsName, skippedNames,
				abundInRow, occuencesInRow,
				NoOfMatrices, writeFiles, true);

	rareStruct* tmpCDR 				= new rareStruct();// 	= {*div, retCnts};
	tmpCDR->div 				= div;
	tmpCDR->cnts 			= cnts;
	tmpCDR->cntsName 	= cntsName;
	tmpCDR->skippedNames		= skippedNames;

	delete cur;
	return tmpCDR;
}


rareStruct* calcDivRarVec(int i, vector<string> fileNames, DivEsts* div, long rareDep,
	vector<vector<uint>>* abundInRow, vector<vector<uint>>* occuencesInRow,
	string outF, int repeats, int writeFiles){
//	Rcout << i << " ";
	smplVec* cur = new smplVec(fileNames[i],4);

	//div->SampleName = curS;
	std::vector<vector<uint>> cnts;
	vector<rare_map> cntsMap;
	string cntsName;
	string skippedNames;
	cur->rarefy(rareDep, outF, repeats,
					div, cntsMap, cntsName, skippedNames, abundInRow, occuencesInRow,
					writeFiles, false,writeFiles);

	//delete cur;
	//return div;
	rareStruct* tmpRS 			= new rareStruct();// 	= {*div, retCnts};
	tmpRS->div 							= div;
	tmpRS->cnts 						= cntsMap;
	tmpRS->cntsName 				= cntsName;
	tmpRS->skippedNames			= skippedNames;

	delete cur;

	//if( remove( fileNames[i].c_str() ) != 0 ){
		//cerr << "LowMem: Error deleting file: " << fileNames[i] << std::endl;
	//}
	return tmpRS;
}


void helpMsg(){
}



void rareLowMem(string inF, string outF, int NoOfMatrices, long arg4, int repeats,
	vector<DivEsts*> *  divvs,
	std::vector<vector<rare_map>>& MaRare,
	std::vector<string>& cntsNames,
	std::vector<string>& skippedNames,
	std::vector<string>& rowNames, int numThr, bool verbose ){
	// this mode takes the file, reads it in memory
	// prints the columns to their own files
	// then it loads those files again and
	// rarefies each column
	// the measures are then combines again.

	//split mat code

	vector<string> fileNames;
	Matrix* Mo 	= new Matrix(inF, outF, "", fileNames, false, true);
	vector < string > SampleNames 	= Mo->getSampleNames();
	rowNames 		= Mo->getRowNames();

	int rareDep 	= arg4;
	if(rareDep == 0){
		// rarefy to smallest colSum
		rareDep = round(0.95 * Mo->getMinColSum());
		//if(rareDep == 0){
		//	cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
		//	return;
		//}
	}
	delete Mo;

	int done = 0; // number of samples processed for multithreading
	uint i = 0;
	std::future<rareStruct*> *tt = new std::future<rareStruct*>[numThr - 1];
	// abundance vectors to hold the number of occurences of genes per row
	// this will be used for Chao2 estimation
	vector<vector<uint>> abundInRow(repeats, vector<uint>(Mo->rowNum(),0));
	vector<vector<uint>> occuencesInRow(repeats, vector<uint>(Mo->rowNum(),0));

	//rarefection code
	while(i < fileNames.size()){

		// allow multithreading
		//if(verbose == true){
		//	cerr << "At Sample " << i+1 << " of " << fileNames.size() << " Samples";
		//}
		uint toWhere = done + numThr - 1;
		if ((uint)((uint)fileNames.size() - 2 ) < toWhere){
			toWhere = fileNames.size() - 2;
		}
		// launch samples in threads
		for (; i < toWhere; i++){
			DivEsts * div 	= new DivEsts();
			div->SampleName = SampleNames[i];
			tt[i - done] = async(std::launch::async, calcDivRarVec, i, fileNames,  div, rareDep,
				 									&abundInRow, &occuencesInRow, outF, repeats, NoOfMatrices);
		}
		// launch one in the mainthread
		DivEsts * div 	= new DivEsts();
		div->SampleName = SampleNames[i];
		rareStruct* tmpRS;
		tmpRS = calcDivRarVec(i, fileNames,  div, rareDep,
													 &abundInRow, &occuencesInRow, outF, repeats, NoOfMatrices);
		i++;

		// process created data, first threads, then main thread
		i = done;
		for (; i < toWhere; i++){
			rareStruct* RSasync;
			RSasync 		= tt[i-done].get();
			divvs->push_back(RSasync->div);
			string curS 	= SampleNames[i];
			//divvs[i-done]->print2file(outF + curS + "_alpha_div.tsv");

			// add the matrices to the container
			if(NoOfMatrices > 0){
				for(uint i = 0; i < RSasync->cnts.size(); i++){
					MaRare[i].push_back(RSasync->cnts[i]);
				}
				// save sample name for naming purposes
				if(RSasync->cntsName.size() != 0){
					cntsNames.push_back(RSasync->cntsName);
				}
			}
			delete RSasync;
		}

		// main thread divv push back
		divvs->push_back(tmpRS->div);
		string curS 	= SampleNames[i];
		//divvs[i]->print2file(outF + curS + "_alpha_div.tsv");
		if(NoOfMatrices > 0){
			for(uint i = 0; i < tmpRS->cnts.size(); i++){
				MaRare[i].push_back(tmpRS->cnts[i]);
			}

			// save sample name for naming purposes
			if(tmpRS->cntsName.size() != 0){
				cntsNames.push_back(tmpRS->cntsName);
			}
		}
		delete tmpRS;
		i++;
		done = i;
	}

}

















//int main(int argc, char* argv[])
int rarefyMain(string inF, string outF, string mode,
	int repeats, long rareDep, unsigned int numThr , bool verbose,
	vector<vector<mat_fl>> rmatrix,
	vector< string > cnames , vector< string > rnames ,
	vector<DivEsts*> *  divvs,
	std::vector<vector<rare_map>> &retCnts,
	std::vector<string>& cntsNames,
	std::vector<string>& skippedNames,
	std::vector<mat_fl>& ACE,
	std::vector<mat_fl>& ICE,
	std::vector<mat_fl>& chao2,
	std::vector<string>& rowNames, int NoOfMatrices,
	bool transpose)
{
	// compatibility to main rare software
	bool writeFiles = false;

	MyRNG rng;


	if (mode == "rare_inmat"){
		Matrix* Mo;
		if(inF != ""){
			  Mo = new Matrix(inF, "");//no arg for outfile &  hierachy | gene subset
		}else{
		  Mo = new Matrix(); // empty matrix for filling later
		  int nr = rmatrix.size();
		  for(int i=0; i < nr; i++){
		    Mo->addRow(rmatrix[i]);
		  }
		  Mo->setSampleNames(cnames);
		  Mo->setRowNames(rnames);
		}

		if(transpose == true){

			Mo->transpose();

		}

		if(rareDep == 0){
			// rarefy to smallest colSum
			rareDep = round(0.95 * Mo->getMinColSum());
			//if(rareDep == 0){
			//	cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
			//	return 0;
			//}
		}
		rowNames = Mo->getRowNames();
		// abundance vectors to hold the number of occurences of genes per row
		// this will be used for Chao2 estimation
		vector<vector<uint>> abundInRow(repeats, vector<uint>(Mo->rowNum(),0));
		vector<vector<uint>> occuencesInRow(repeats, vector<uint>(Mo->rowNum(),0));

		//vector<DivEsts*> divvs(Mo->smplNum(),NULL);
		//divvs->resize(0,NULL); // paul resize the vector

		//cerr << "TH";
		//std::future<DivEsts*> *tt = new std::future<DivEsts*>[numThr - 1];
		std::future<rareStruct*> *tt = new std::future<rareStruct*>[numThr - 1];
		//Rcout << "threads\n";
		uint i = 0; uint done = 0;
		while ( i < Mo->smplNum()){

			uint toWhere = done+numThr - 1;
			if ((uint)((uint)Mo->smplNum() - 2 ) < toWhere){
					toWhere = Mo->smplNum() - 2;
			}
			for (; i < toWhere; i++){ // with just one thread this is not used?
				DivEsts * div = new DivEsts();
				tt[i - done] = async(std::launch::async, calcDivRar, i, Mo, div, rareDep,
								 	&abundInRow, &occuencesInRow, outF,
									repeats, writeFiles, NoOfMatrices);

			}

			//use main thread to calc one sample as well
			DivEsts * div = new DivEsts();
			//divvs[i] = calcDivRar(i, Mo, div, rareDep, outF, repeats, writeFiles);
			rareStruct* tmpCDr;
			//tmpCDr = new rareStruct;
			tmpCDr 		= calcDivRar(i, Mo, div, rareDep,  &abundInRow, &occuencesInRow,
									outF, repeats, writeFiles,
									 NoOfMatrices);



			//Rcout <<'\n' <<  i << " retCnts: " << tmpCDr->retCnts.size();
			i++;
			i 			= done;
			for (; i < toWhere; i++){
				rareStruct* CDrAsync;
				CDrAsync = tt[i-done].get();
				// append diversity measures
				divvs->push_back(CDrAsync->div);

				if(NoOfMatrices > 0){
					// append vector to matrix
					uint repI = 0;
					while(repI < CDrAsync->cnts.size()){
						retCnts[repI].push_back(CDrAsync->cnts[repI]);
						repI++;
					}
					// save sample name for naming purposes
					if(CDrAsync->cntsName.size() != 0){
						cntsNames.push_back(CDrAsync->cntsName);
					}
				}
				// skippedNames
				if(CDrAsync->skippedNames.size() > 0){
					skippedNames.push_back(CDrAsync->skippedNames);
				}
				delete CDrAsync;
			}

			// main thread
			divvs->push_back(tmpCDr->div);
			if(NoOfMatrices > 0){
				uint repI = 0;
				while(repI < tmpCDr->cnts.size()){
					retCnts[repI].push_back(tmpCDr->cnts[repI]);
					repI++;
				}
				// save sample name for naming purposes
				if(tmpCDr->cntsName.size() != 0){
					cntsNames.push_back(tmpCDr->cntsName);
				}
			}
			// skippedNames
			if(tmpCDr->skippedNames.size() > 0){
				skippedNames.push_back(tmpCDr->skippedNames);
			}


			delete tmpCDr;



			i++;
			done = i;
		}
		delete Mo;

		// compute chao2, ACE, ICE and write to file
		computeChao2(chao2, abundInRow);
		computeCE(ICE, abundInRow);
		computeCE(ACE, occuencesInRow);


	}else if(mode == "rare_lowMem"){
		rareLowMem(inF, outF, NoOfMatrices,  rareDep,  repeats,
		divvs, retCnts, cntsNames, skippedNames, rowNames, numThr, verbose);
	}


	return 0;
}
