
#include "RRare.h"


const char* rar_ver="0.63		 alpha R";

struct cDR{
	DivEsts* div;
	std::vector<vector<vector<uint>>> retCnts;
	string retCntsSampleName;
	vector<map<uint, uint>> RareSample;
	string skippedSample;
};

cDR* calcDivRar(int i, Matrix* Mo, DivEsts* div,  long rareDep, string outF,
	int repeats, int writeFiles,
	int NoOfMatrices){

	smplVec* cur 	= Mo->getSampleVec(i);
	string curS 	= Mo->getSampleName(i);
	div->SampleName = curS;

	// vector holding the rarefaction results for this sample
	// repeat times
	std::vector<map<uint, uint>> RareSample;
	string retCntsSampleName;
	string skippedSample;
	cur->rarefy(rareDep, outF, repeats, div, RareSample, retCntsSampleName, skippedSample,
				NoOfMatrices, writeFiles, true);

	cDR* tmpCDR 				= new cDR();// 	= {*div, retCnts};
	tmpCDR->div 				= div;
	tmpCDR->RareSample 			= RareSample;
	tmpCDR->retCntsSampleName 	= retCntsSampleName;
	tmpCDR->skippedSample		= skippedSample;

	delete cur;
	return tmpCDR;
}


void helpMsg(){
	string  AvailableModes = "Available run modes:\nnormalize\nsplitMat\nlineExtr\nmergeMat\nsumMat\nrarefaction\nrare_inmat\nmodule\n";
	cerr << AvailableModes << "Provide two arguments\nexiting..\n";
	cerr << "------------------------------\nAuthor: falk.hildebrand@gmail.com\n";
	std::exit(2);
}



void rareLowMem(string inF, string outF, int writeFiles, long arg4, int repeats){
	// this mode takes the file, reads it in memory
	// prints the columns to their own files
	// then it loads those files again and
	// rarefies each column
	// the measures are then combines again.

	//split mat code
	vector<string> fileNames;
	Matrix* Mo 	= new Matrix(inF, outF, "", fileNames, false, true);
	vector < string > SampleNames 	= Mo->getSampleNames();
	vector < string > rowNames 		= Mo->getRowNames();

	int rareDep 	= arg4;
	if(rareDep == 0){
		// rarefy to smallest colSum
		rareDep = round(0.95 * Mo->getMinColSum());
	}
	delete Mo;


	int NoOfMatrices = writeFiles;
	vector< vector< vector< uint > > > MaRare (NoOfMatrices);
	std::vector<string> cntsNames;

	//rarefection code
	vector<DivEsts*> divvs(fileNames.size(),NULL);
	for(uint i = 0; i < fileNames.size(); i++){
		smplVec* cur 		= new smplVec(fileNames[i], 4);
		DivEsts * div 		= new DivEsts();
		div->SampleName 	= SampleNames[i];
		//placeholder for R function, not to be filled here
		std::vector<map<uint, uint>> cnts;
		string cntsName;
		string skippedSample;
		cur->rarefy(rareDep,outF,repeats,div, cnts, cntsName, skippedSample, writeFiles,false,NoOfMatrices);
		divvs[i] 			= div;



		if(NoOfMatrices > 0){
			vector < string > rowIDs = cur->getRowNames();
			vector < uint > nrowIDs(rowIDs.size());
			// convert ids into integer vector
			for(uint i = 0; i < rowIDs.size(); i++){
				nrowIDs[i] = std::stoi(rowIDs[i]);
			}
			for(uint i = 0; i < cnts.size(); i++){
				// reshape each vector, as some are zero, and we need to rematch values and rows
				vector <uint> tmpVec(rowNames.size(), 0);
				for(uint j = 0; j < nrowIDs.size(); j++){
					tmpVec[nrowIDs[j]] = cnts[i][j];
				}
				MaRare[i].push_back(tmpVec);
			}
			// save sample name for naming purposes
			if(cntsName.size() != 0){
				cntsNames.push_back(cntsName);
			}
		}

		delete cur;
	}
	// print the div estimates out into a file
//	printDivMat(outF + "all.txt", divvs);
	for (size_t i = 0; i < divvs.size(); i++){
		delete divvs[i];
	}
	if(NoOfMatrices > 0){
		for(uint i = 0; i < MaRare.size(); i++){
		//	printRareMat(outF + "rarefied_" +  std::to_string(i) + ".tsv", MaRare[i], cntsNames, rowNames);
		}
	}

	// delete tmp file we created
	for(uint i = 0; i < fileNames.size(); i++){
		if( remove( fileNames[i].c_str() ) != 0 ){
			cerr << "Error deleting file: " << fileNames[i];
		}
	}

	cout << "Finished\n";
}

















//int main(int argc, char* argv[])
int rarefyMain(string inF, string mode,
	int repeats, long rareDep, unsigned int numThr , bool verbose,
	vector<vector<mat_fl>> rmatrix,
	vector< string > cnames , vector< string > rnames ,
	vector<DivEsts*> *  divvs,
	std::vector<vector<map<uint, uint>>> &retCnts,
	std::vector<string>& retCntsSampleNames,
	std::vector<string>& skippedSamples,
	std::vector<string>& rowNames, int NoOfMatrices,
	bool transpose)
{
	// compatibility to main rare software
	bool writeFiles = false;
	string outF = "";

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
		}
		rowNames = Mo->getRowNames();
		//vector<DivEsts*> divvs(Mo->smplNum(),NULL);
		//divvs->resize(0,NULL); // paul resize the vector
		if(verbose == true){
			cout << "Using " << numThr << " threads\n";
		}
		//cerr << "TH";
		//std::future<DivEsts*> *tt = new std::future<DivEsts*>[numThr - 1];
		std::future<cDR*> *tt = new std::future<cDR*>[numThr - 1];
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
			cDR* tmpCDr;
			//tmpCDr = new cDR;
			tmpCDr 		= calcDivRar(i, Mo, div, rareDep, outF, repeats, writeFiles,
									 NoOfMatrices);



			//cout <<'\n' <<  i << " retCnts: " << tmpCDr->retCnts.size();
			i++;
			i 			= done;
			for (; i < toWhere; i++){
				cDR* CDrAsync;
				CDrAsync = tt[i-done].get();
				// append diversity measures
				divvs->push_back(CDrAsync->div);

				if(NoOfMatrices > 0){
					// append vector to matrix
					int repI = 0;
					while(repI < CDrAsync->RareSample.size()){
						retCnts[repI].push_back(CDrAsync->RareSample[repI]);
						repI++;
					}
					// save sample name for naming purposes
					if(CDrAsync->retCntsSampleName.size() != 0){
						retCntsSampleNames.push_back(CDrAsync->retCntsSampleName);
					}
				}
				delete CDrAsync;
			}

			// main thread
			divvs->push_back(tmpCDr->div);
			if(NoOfMatrices > 0){
				int repI = 0;
				while(repI < tmpCDr->RareSample.size()){
					retCnts[repI].push_back(tmpCDr->RareSample[repI]);
					repI++;
				}
				// save sample name for naming purposes
				if(tmpCDr->retCntsSampleName.size() != 0){
					retCntsSampleNames.push_back(tmpCDr->retCntsSampleName);
				}



			}


			delete tmpCDr;



			i++;
			done = i;
		}
		delete Mo;
	}else if(mode == "rare_lowMem"){
		rareLowMem(inF, outF, writeFiles,  rareDep,  repeats);
	}


	return 0;
}
