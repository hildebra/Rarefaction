
#include "RRare.h"


const char* rar_ver="0.63		 alpha R";

struct cDR{
	DivEsts* div;
	std::vector<vector<vector<uint>>> retCnts;
	string retCntsSampleName;
	vector<map<uint, uint>> RareSample;
	string skippedSample;
};

struct rareStruct{
	DivEsts* div;
	string cntsName;
	vector< map< uint, uint>> cnts;
	string skippedNames;
};

rareStruct* calcDivEst(int i, Matrix* Mo, DivEsts* div, long rareDep, string outF, int writeFiles){
	smplVec* cur = Mo->getSampleVec(i);
	string curS = Mo->getSampleName(i);
	div->SampleName = curS;
	std::vector<map<uint, uint>> cnts;
	string cntsName;
	string skippedNames;
	vector<mat_fl> flcnts = Mo->getColumn(i);
	cur->estimateDiv(rareDep, div, flcnts,  cnts, cntsName, skippedNames);
	//delete cur;
	//return div;
	rareStruct* tmpRS 			= new rareStruct();// 	= {*div, retCnts};
	tmpRS->div 							= div;
	tmpRS->cnts 						= cnts;
	tmpRS->cntsName 				= cntsName;
	tmpRS->skippedNames			= skippedNames;
	delete cur;
	return tmpRS;
}

rareStruct* calcDivRar(int i, Matrix* Mo, DivEsts* div,  long rareDep, string outF,
	int repeats, int writeFiles,
	int NoOfMatrices){

	smplVec* cur 	= Mo->getSampleVec(i);
	string curS 	= Mo->getSampleName(i);
	div->SampleName = curS;

	// vector holding the rarefaction results for this sample
	// repeat times
	std::vector<map<uint, uint>> cnts;
	string cntsName;
	string skippedNames;
	cur->rarefy(rareDep, outF, repeats, div, cnts, cntsName, skippedNames,
				NoOfMatrices, writeFiles, true);

	rareStruct* tmpCDR 				= new rareStruct();// 	= {*div, retCnts};
	tmpCDR->div 				= div;
	tmpCDR->cnts 			= cnts;
	tmpCDR->cntsName 	= cntsName;
	tmpCDR->skippedNames		= skippedNames;

	delete cur;
	return tmpCDR;
}


rareStruct* calcDivRarVec(int i, vector<string> fileNames, DivEsts* div, long rareDep, string outF, int repeats, int writeFiles){
//	cout << i << " ";
	smplVec* cur = new smplVec(fileNames[i],4);

	//div->SampleName = curS;
	std::vector<vector<uint>> cnts;
	vector< map< uint, uint>> cntsMap;
	string cntsName;
	string skippedNames;
	cur->rarefy(rareDep, outF, repeats,
					div, cntsMap, cntsName, skippedNames,
					writeFiles, false,writeFiles);

	//delete cur;
	//return div;
	rareStruct* tmpRS 			= new rareStruct();// 	= {*div, retCnts};
	tmpRS->div 							= div;
	tmpRS->cnts 						= cntsMap;
	tmpRS->cntsName 				= cntsName;
	tmpRS->skippedNames			= skippedNames;

	delete cur;

	if( remove( fileNames[i].c_str() ) != 0 ){
		cerr << "LowMem: Error deleting file: " << fileNames[i] << std::endl;
	}
	return tmpRS;
}


void helpMsg(){
	string  AvailableModes = "Available run modes:\nnormalize\nsplitMat\nlineExtr\nmergeMat\nsumMat\nrarefaction\nrare_inmat\nmodule\n";
	cerr << AvailableModes << "Provide two arguments\nexiting..\n";
	cerr << "------------------------------\nAuthor: falk.hildebrand@gmail.com\n";
	std::exit(2);
}



void rareLowMem(string inF, string outF, int NoOfMatrices, long arg4, int repeats,
	vector<DivEsts*> *  divvs,
	std::vector<vector<map<uint, uint>>>& MaRare,
	std::vector<string>& cntsNames,
	std::vector<string>& skippedNames,
	std::vector<string>& rowNames, int numThr, bool verbose ){
	// this mode takes the file, reads it in memory
	// prints the columns to their own files
	// then it loads those files again and
	// rarefies each column
	// the measures are then combines again.

	//split mat code
	if(verbose == true){
		cout << " Low mem test" << std::endl;
		cout << "Tmp dir: " << outF << std::endl;
	}
	vector<string> fileNames;
	Matrix* Mo 	= new Matrix(inF, outF, "", fileNames, false, true);
	vector < string > SampleNames 	= Mo->getSampleNames();
	rowNames 		= Mo->getRowNames();

	int rareDep 	= arg4;
	if(rareDep == 0){
		// rarefy to smallest colSum
		rareDep = round(0.95 * Mo->getMinColSum());
		if(rareDep == 0){
			cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
			return;
		}
	}
	delete Mo;

	int done = 0; // number of samples processed for multithreading
	uint i = 0;
	std::future<rareStruct*> *tt = new std::future<rareStruct*>[numThr - 1];

	//int NoOfMatrices = writeFiles;

	//rarefection code
	//divvs->resize(fileNames.size());
	while(i < fileNames.size()){

		// allow multithreading
		if(verbose == true){
			cerr << "At Sample " << i+1 << " of " << fileNames.size() << " Samples";
		}
		uint toWhere = done + numThr - 1;
		if ((uint)((uint)fileNames.size() - 2 ) < toWhere){
			toWhere = fileNames.size() - 2;
		}
		// launch samples in threads
		for (; i < toWhere; i++){
			DivEsts * div 	= new DivEsts();
			div->SampleName = SampleNames[i];
			tt[i - done] = async(std::launch::async, calcDivRarVec, i, fileNames,  div, rareDep, outF, repeats, NoOfMatrices);
		}
		// launch one in the mainthread
		DivEsts * div 	= new DivEsts();
		div->SampleName = SampleNames[i];
		rareStruct* tmpRS;
		tmpRS = calcDivRarVec(i, fileNames,  div, rareDep, outF, repeats, NoOfMatrices);
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


	// delete tmp file we created
	fileNames.push_back(outF + "sums.txt");
	/*for(uint i = 0; i < fileNames.size(); i++){
		if( remove( fileNames[i].c_str() ) != 0 ){
			cerr << "Error deleting file: " << fileNames[i];
		}
	}*/

}





int estimateMode(string inF, string outF, int NoOfMatrices,
									uint rareDep, vector<DivEsts*> *  divvs ,
									std::vector<vector<map<uint, uint>>> &retCnts,
									vector<vector<mat_fl>>  rmatrix,
									std::vector<string>& cntsNames,
									std::vector<string>& skippedNames,
									vector< string > rnames,
									vector< string > cnames,
									bool transpose,
									uint numThr, bool verbose){

	int writeFiles = 0;


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
		if(verbose == true){
			cout << "Will now transpose the matrix\n";
		}
		Mo->transpose();
		if(verbose == true){
			cout << "Done transposing\n";
		}
	}

	if(rareDep == 0){
		// rarefy to smallest colSum
		rareDep = round(0.95 * Mo->getMinColSum());
		if(rareDep == 0){
			cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
			return 0;
		}
	}

	//Matrix* Mo = new Matrix(inF, "");

	//vector<DivEsts*> divvs(Mo->smplNum(),NULL);

	Mo->normalize(rareDep, false);

	uint i = 0; uint done = 0;
	while ( i < Mo->smplNum()){

		rareStruct* tmpRS;
		DivEsts * div 	= new DivEsts();
		tmpRS = calcDivEst(i, Mo, div, rareDep, "", writeFiles);

		divvs->push_back(tmpRS->div);

		if(NoOfMatrices > 0){
			uint repI = 0;
			while(repI < tmpRS->cnts.size()){
				cout << "tmpRs size " << tmpRS->cnts.size() << "repI  " <<repI << std::endl;
				retCnts[repI].push_back(tmpRS->cnts[repI]);
				repI++;
			}
			// save sample name for naming purposes
			if(tmpRS->cntsName.size() != 0){
				cntsNames.push_back(tmpRS->cntsName);
			}
		}
		// skippedNames
		if(tmpRS->skippedNames.size() > 0){
			skippedNames.push_back(tmpRS->skippedNames);
		}

		i++;
	}
	//Mo->estimateDiversity(Depth);

	//printDivMat(outF , divvs, true);
	//Mo->writeMatrix(outF);
	delete Mo;
}













//int main(int argc, char* argv[])
int rarefyMain(string inF, string outF, string mode,
	int repeats, long rareDep, unsigned int numThr , bool verbose,
	vector<vector<mat_fl>> rmatrix,
	vector< string > cnames , vector< string > rnames ,
	vector<DivEsts*> *  divvs,
	std::vector<vector<map<uint, uint>>> &retCnts,
	std::vector<string>& cntsNames,
	std::vector<string>& skippedNames,
	std::vector<string>& rowNames, int NoOfMatrices,
	bool transpose)
{
	// compatibility to main rare software
	bool writeFiles = false;

	MyRNG rng;


	if (mode == "memory"){
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
			if(verbose == true){
				cout << "Will now transpose the matrix\n";
			}
			Mo->transpose();
			if(verbose == true){
				cout << "Done transposing\n";
			}
		}

		if(rareDep == 0){
			// rarefy to smallest colSum
			rareDep = round(0.95 * Mo->getMinColSum());
			if(rareDep == 0){
				cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
				return 0;
			}
		}
		rowNames = Mo->getRowNames();
		//vector<DivEsts*> divvs(Mo->smplNum(),NULL);
		//divvs->resize(0,NULL); // paul resize the vector
		if(verbose == true){
			cout << "Using " << numThr << " threads\n";
		}
		//cerr << "TH";
		//std::future<DivEsts*> *tt = new std::future<DivEsts*>[numThr - 1];
		std::future<rareStruct*> *tt = new std::future<rareStruct*>[numThr - 1];
		//cout << "threads\n";
		uint i = 0; uint done = 0;
		while ( i < Mo->smplNum()){
			if(verbose == true){
				cout << "At Sample " << i+1 << " of " << Mo->smplNum() << " Samples\n";
			}
			uint toWhere = done+numThr - 1;
			if ((uint)((uint)Mo->smplNum() - 2 ) < toWhere){
					toWhere = Mo->smplNum() - 2;
			}
			for (; i < toWhere; i++){ // with just one thread this is not used?
				DivEsts * div = new DivEsts();
				tt[i - done] = async(std::launch::async, calcDivRar, i, Mo, div, rareDep, outF,
									repeats, writeFiles, NoOfMatrices);

			}

			//use main thread to calc one sample as well
			DivEsts * div = new DivEsts();
			//divvs[i] = calcDivRar(i, Mo, div, rareDep, outF, repeats, writeFiles);
			rareStruct* tmpCDr;
			//tmpCDr = new rareStruct;
			tmpCDr 		= calcDivRar(i, Mo, div, rareDep, outF, repeats, writeFiles,
									 NoOfMatrices);



			//cout <<'\n' <<  i << " retCnts: " << tmpCDr->retCnts.size();
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
	}else if(mode == "swap"){
		rareLowMem(inF, outF, NoOfMatrices,  rareDep,  repeats,
		divvs, retCnts, cntsNames, skippedNames, rowNames, numThr, verbose);
	}else if(mode == "estimate"){


		estimateMode(inF, outF, NoOfMatrices, rareDep, divvs ,
			retCnts,  rmatrix,cntsNames,
			skippedNames, rnames, cnames,
			transpose, numThr,  verbose);


		//estimateMode(inF, outF, rareDep, divvs, retCnts, );

	}


	return 0;
}
