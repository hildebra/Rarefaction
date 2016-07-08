// rarefaction.cpp
// usage: rarefaction in_matrix outfile
//outfile contains richness and sample sum; in same dir several files are created
//C:\Users\Falk\SkyDrive\science\data\test\test.mat.subs

//#include "Matrix.h"
#include "ClStr2Mat.h"
#include "Rare.h"

const char* rar_ver="0.64 alpha";


rareStruct* calcDivRar(int i, Matrix* Mo, DivEsts* div, long rareDep, string outF, int repeats, int writeFiles){
	smplVec* cur = Mo->getSampleVec(i);
	string curS = Mo->getSampleName(i);
	div->SampleName = curS;
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
	return tmpRS;
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
	tmpRS->IDs 							= cur->getRowNames();

	delete cur;
	return tmpRS;
}



void helpMsglegacy(){
	string  AvailableModes = "Available run modes:\nnormalize\nsplitMat\nlineExtr\nmergeMat\nsumMat\nrarefaction\nrare_inmat\nmodule\n";
	cerr << AvailableModes << "Provide two arguments\nexiting..\n";
	cerr << "------------------------------\nAuthor: falk.hildebrand@gmail.com\n";
	std::exit(2);
}

void stateVersion(){
	printf("rare %s\n", rar_ver);
}

void helpMsg(){
	stateVersion();
	printf("\n");
	printf("usage: rare mode input output [options] \n");
	//printf("\n");
	printf("available run modes:\n");
	//printf("    normalize     TODO\n");
	//printf("    splitMat      TODO\n");
	//printf("    lineExtr      TODO\n");
	//printf("    mergeMat      TODO\n");
	//printf("    sumMat        TODO\n");
	printf("    rarefaction   TODO\n");
	printf("    rare_lowMem   rarefy a matrix and compute diversity measures. Using less memory.\n");
	printf("    rare_inmat    rarefy a matrix and compute diversity measures. Process in memory.\n");
	//printf("    module        TODO\n");
	printf("\n");
	printf("details:\n");
	printf("\n");
	printf("\n");
	printf("mode:  rare_lowMem\n");
	printf("usage: rare rare_lowMem input output depth repeats write threads\n");
	printf("\n");
	printf("       mode      rare_inmat (Other modes are not available at the moment. Can be any string.)\n");
	printf("       input     path to a .csv file                   (required)\n");
	printf("       output    path for the ouput and tmp files      (required)\n");
	printf("       depth     rarefaction depth                     (default: 1000)\n");
	printf("       repeats   number of times to compute diversity  (default: 1)\n");
	printf("       write     number of files to write              (default: same as repeats)\n");
	printf("       threads   number of threads to use in parallel  (default: 1)\n");
	printf("\n");
	printf("mode:  rare_inmat\n");
	printf("usage: rare rare_inmat input output depth repeats write threads\n");
	printf("\n");
	printf("       mode      rare_inmat (Other modes are not available at the moment. Can be any string.)\n");
	printf("       input     path to a .csv file                   (required)\n");
	printf("       output    path for the ouput                    (required)\n");
	printf("       depth     rarefaction depth                     (default: 1000)\n");
	printf("       repeats   number of times to compute diversity  (default: 1)\n");
	printf("       write     number of files to write              (default: same as repeats)\n");
	printf("       threads   number of threads to use in parallel  (default: 1)\n");
	printf("\n");


	std::exit(2);
}

/*

void rareLowMem(string inF, string outF, int writeFiles, string arg4, int repeats, int numThr = 1){
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

	int rareDep 	= atoi(arg4.c_str());
	if(rareDep == 0){
		// rarefy to smallest colSum
		rareDep = round(0.95 * Mo->getMinColSum());
		if(rareDep == 0.0){
			cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
			exit(1);
		}
	}
	delete Mo;


	int NoOfMatrices = writeFiles;
	vector< vector< map< uint, uint > > > MaRare (NoOfMatrices);
	std::vector<string> cntsNames;

	int done = 0; // number of samples processed for multithreading
	uint i = 0;
	std::future<rareStruct*> *tt = new std::future<rareStruct*>[numThr - 1];

	//rarefection code
	vector<DivEsts*> divvs(fileNames.size(),NULL);
	while(i < fileNames.size()){

		// allow multithreading
		cerr << "At Sample " << i << " of " << fileNames.size() << " Samples";
		uint toWhere = done + numThr - 1;
		if ((uint)((uint)fileNames.size() - 2 ) < toWhere){
			toWhere = fileNames.size() - 2;
		}
		// launch samples in threads
		for (; i < toWhere; i++){
			DivEsts * div 	= new DivEsts();
			div->SampleName = SampleNames[i];
			tt[i - done] = async(std::launch::async, calcDivRarVec, i, fileNames,  div, rareDep, outF, repeats, writeFiles);
		}
		// launch one in the mainthread
		DivEsts * div 	= new DivEsts();
		div->SampleName = SampleNames[i];
		rareStruct* tmpRS;
		tmpRS = calcDivRarVec(i, fileNames,  div, rareDep, outF, repeats, writeFiles);
		i++;

		// process created data, first threads, then main thread
		i = done;
		for (; i < toWhere; i++){
			rareStruct* RSasync;
			RSasync 		= tt[i-done].get();
			divvs[i] 		= RSasync->div;
			string curS 	= SampleNames[i];
			divvs[i-done]->print2file(outF + curS + "_alpha_div.tsv");

			// add the matrices to the container
			if(NoOfMatrices > 0){
				vector < string > rowIDs = RSasync->IDs;
				vector < uint > nrowIDs(rowIDs.size());
				// convert ids into integer vector
				for(uint i = 0; i < rowIDs.size(); i++){
					nrowIDs[i] = std::stoi(rowIDs[i]);
				}
				for(uint i = 0; i < RSasync->cnts.size(); i++){
					// reshape each vector, as some are zero, and we need to rematch values and rows
					// we use the row Ids which we created correctly when splitting the vector from the input file
					std::map <uint, uint> tmpVec;
						for (auto const& x : RSasync->cnts[i]){
							tmpVec[nrowIDs[x.first]] = x.second;
						}
					MaRare[i].push_back(tmpVec);
				}
				//for(uint i = 0; i < RSasync->cnts.size(); i++){
				//	MaRare[i].push_back(RSasync->cnts[i]);
				//}
				// save sample name for naming purposes
				if(RSasync->cntsName.size() != 0){
					cntsNames.push_back(RSasync->cntsName);
				}
			}
			delete RSasync;
		}

		// main thread divv push back
		divvs[i] = tmpRS->div;
		string curS 	= SampleNames[i];
		divvs[i]->print2file(outF + curS + "_alpha_div.tsv");

		if(NoOfMatrices > 0){
					vector < string > rowIDs = tmpRS->IDs;
					vector < uint > nrowIDs(rowIDs.size());
					// convert ids into integer vector
					for(uint i = 0; i < rowIDs.size(); i++){
						nrowIDs[i] = std::stoi(rowIDs[i]);
					}
					for(uint i = 0; i < tmpRS->cnts.size(); i++){
						// reshape each vector, as some are zero, and we need to rematch values and rows
						// we use the row Ids which we created correctly when splitting the vector from the input file
						std::map <uint, uint> tmpVec;
							for (auto const& x : tmpRS->cnts[i]){
								tmpVec[nrowIDs[x.first]] = x.second;
							}
						MaRare[i].push_back(tmpVec);
					}
					// save sample name for naming purposes
					if(tmpRS->cntsName.size() != 0){
						cntsNames.push_back(tmpRS->cntsName);
					}
				}


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



	// print the div estimates out into a file
	printDivMat(outF + "median_alpha_diversity.tsv", divvs);
	for (size_t i = 0; i < divvs.size(); i++){
		delete divvs[i];
	}
	if(NoOfMatrices > 0){
		for(uint i = 0; i < MaRare.size(); i++){
			printRareMat(outF + "rarefied_to_" + std::to_string(rareDep) + "_n_" +  std::to_string(i) + ".tsv", MaRare[i], cntsNames, rowNames);
		}
	}

	// delete tmp file we created
	for(uint i = 0; i < fileNames.size(); i++){
		if( remove( fileNames[i].c_str() ) != 0 ){
			cerr << "Error deleting file: " << fileNames[i] << std::endl;
		}
	}

	cout << "Finished\n";


}

*/

void rareExtremLowMem(string inF, string outF, int writeFiles, string arg4, int repeats, int numThr = 1, bool storeBinary = false){
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

	int rareDep 	= atoi(arg4.c_str());
	if(rareDep == 0){
		// rarefy to smallest colSum
		rareDep = round(0.95 * Mo->getMinColSum());
		if(rareDep == 0.0){
			cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
			exit(1);
		}
	}
	delete Mo;


	int NoOfMatrices = writeFiles;
	vector< vector< map< uint, uint > > > MaRare (NoOfMatrices);
	std::vector<string> cntsNames;
	vector < vector < string > > tmpMatFiles (NoOfMatrices );
	int done = 0; // number of samples processed for multithreading
	uint i = 0;
	std::future<rareStruct*> *tt = new std::future<rareStruct*>[numThr - 1];


	//rarefection code
	vector<DivEsts*> divvs(fileNames.size(),NULL);
	while(i < fileNames.size()){

		// allow multithreading
		cerr << "At Sample " << i << " of " << fileNames.size() << " Samples";
		uint toWhere = done + numThr - 1;
		if ((uint)((uint)fileNames.size() - 2 ) < toWhere){
			toWhere = fileNames.size() - 2;
		}
		// launch samples in threads
		for (; i < toWhere; i++){
			DivEsts * div 	= new DivEsts();
			div->SampleName = SampleNames[i];
			tt[i - done] = async(std::launch::async, calcDivRarVec, i, fileNames,  div, rareDep, outF, repeats, writeFiles);
		}
		// launch one in the mainthread
		DivEsts * div 	= new DivEsts();
		div->SampleName = SampleNames[i];
		rareStruct* tmpRS;
		tmpRS = calcDivRarVec(i, fileNames,  div, rareDep, outF, repeats, writeFiles);
		i++;

		// process created data, first threads, then main thread
		i = done;
		for (; i < toWhere; i++){
			rareStruct* tmpRS;
			tmpRS 		= tt[i-done].get();
			divvs[i] 		= tmpRS->div;
			string curS 	= SampleNames[i];
			divvs[i-done]->print2file(outF + curS + "_alpha_div.tsv");

			// add the matrices to the container
			if(NoOfMatrices > 0){
				if(storeBinary == true){
					binaryStoreSample(tmpMatFiles, tmpRS, rowNames,outF, cntsNames, true);
				}else{
					memoryStoreSample(tmpRS, MaRare, cntsNames, true);
				}
			}

			delete tmpRS;
		}

		// main thread divv push back
		divvs[i] 			= tmpRS->div;
		string curS 	= SampleNames[i];
		divvs[i]->print2file(outF + curS + "_alpha_div.tsv");


		// add the matrices to the container
		if(NoOfMatrices > 0){
			if(storeBinary == true){
				binaryStoreSample(tmpMatFiles, tmpRS, rowNames,outF, cntsNames, true);
			}else{
				memoryStoreSample(tmpRS, MaRare, cntsNames, true);
			}
		}

		delete tmpRS;
		i++;
		done = i;
	}

	// delete tmp files we created for rarefaction, as we now do not need them anymore
	for(uint i = 0; i < fileNames.size(); i++){
		if( remove( fileNames[i].c_str() ) != 0 ){
			cerr << "Error deleting file: " << fileNames[i] << std::endl;
		}
	}

	// print the div estimates out into a file
	printDivMat(outF + "median_alpha_diversity.tsv", divvs);
	for (size_t i = 0; i < divvs.size(); i++){
		delete divvs[i];
	}

	// write rarefaction matrices to disk
	if(NoOfMatrices > 0){
		//vector< string > rowNames = Mo->getRowNames();
		if(storeBinary == true){
			printRarefactionMatrix(tmpMatFiles, outF, rareDep, cntsNames, rowNames);
		}else{
			printRarefactionMatrix(MaRare, outF, rareDep, cntsNames, rowNames);
		}
	}

	cout << "Finished\n";
}






int main(int argc, char* argv[])
{

	if (argc < 3) {
		helpMsg();
	}
	stateVersion();

	long rareDep = 1000;	int repeats (1);
	//bool splitMode(false),mergeMode(false),sumUpMode(false);
	string inF = argv[2];
	string outF = argv[3];
	// hardcode mode, as the others are excluded for newRow
	// can be removed, as soon as other functionality is ready
	// to be used. June 2016
	//string mode = argv[1];
	string mode = argv[1];
	string arg4 = "";
	uint numThr = 1; //number of threads to use


	if (argc>=5){
		arg4 = argv[4];
	}

	if (argc>5){
		repeats = atoi(argv[5]);
	}
	int writeFiles = repeats;
	if (argc>6){
		writeFiles = atoi(argv[6]);
	}
	if (argc > 7){
		numThr = atoi(argv[7]);
	}
	if( writeFiles > repeats){
		cerr << "Can not write more files than repeats. Set writefiles to repeats" << std::endl;
		int writeFiles = repeats;
	}

	// flag to store data on disk
	bool storeBinary = true;
	vector < vector < string > > tmpMatFiles (writeFiles );
	if (argc > 8){
		storeBinary = atoi(argv[8]);
	}

	// start the modes
	if (argc>3){
		if (mode == "splitMat") {
			vector<string> empt;
			Matrix* Mo = new Matrix(inF, outF, arg4, empt, false);
			//Mo->splitOnHDD(outF);	//Mo->writeSums(outF);
			delete Mo;
			std::exit(0);
		//}else if (mode == "rare_lowMem") {
		//	rareLowMem(inF, outF, writeFiles,  arg4,  repeats, numThr);
		}else if (mode == "rare_lowMem") {
			rareExtremLowMem(inF, outF, writeFiles,  arg4,  repeats, numThr, storeBinary);
		}	else if (mode == "correl2"){
			//usage: ./rare correl2 [signature matrix] [output matrix] [big gene matrix]
			//reads in signature matrix (e.g. 40 marker genes)
			SigMatrix* Sig = new SigMatrix(inF);
			//Sig->estimateBinModel();
			//readMatrixLinebyLine(arg4,Sig);

		}else if (mode == "module") {
			if (argc < 7) {
				cout << "Usage: ./rare module [KO matrix] [outputfile] [module DB file] [KO redundancy, int] [Pathway completeness, float 0-1] [Enzyme completeness, float 0-1]\n";
				cerr << "Not enough arguments for \"module\" function\n";
				exit(3);
			}
			Matrix* Mo = new Matrix(inF, ""); //needs to be KO file
			cerr << "Estimate mod AB\n";
			Mo->estimateModuleAbund(argv,argc);// arg4, outF); //arg4 needs to be module file, outF makes the (several) matrices
			delete Mo;
			std::exit(0);
		}
		else if (mode == "normalize") {
			if (argc < 4) {
				cout << "Usage: ./rare normalize [in matrix] [outputfile]\n";
				cerr << "Not enough arguments for \"normalize\" function\n";
				exit(3);
			}
			Matrix* Mo = new Matrix(inF, ""); //needs to be KO file
			Mo->normalize();
			Mo->writeMatrix(outF);
			delete Mo;
			std::exit(0);
		} else if (mode == "help" || mode == "-help" || mode == "-h" || mode == "--help"){
			helpMsg();
		} else if (mode == "lineExtr"){
			lineCntOut(inF, outF, arg4);
			std::exit(0);
		} else if (mode == "mergeMat") {
			VecFiles* VFs = new VecFiles(inF, outF, arg4);
			delete VFs;
			std::exit(0);
		} else if (mode == "sumMat") {
			vector<string> empt;
			Matrix* Mo = new Matrix(inF, outF, arg4, empt, true);
			delete Mo;
			std::exit(0);
		} else if (mode == "rarefaction" || mode == "rare_inmat") {
			rareDep = atoi(arg4.c_str());
		}
		else if (mode == "geneMat"){
			cout << "Gene clustering matrix creation\n";
			if (argc < 5) {cerr << "Needs at least 4 arguments\n"; std::exit(0);}
			ClStr2Mat* cl = new ClStr2Mat(inF,outF,arg4,argv[5]);
			delete cl;
			std::exit(0);
		} else {
			helpMsg();
		}
	}




	MyRNG rng;
	//test rand numbers.. check
	//std::uniform_int_distribution<unsigned long> uint_distx(0,2223951715);	for (int i=0;i<100;i++){cout<<uint_distx(rng)<<" ";}	int i;	cin>>i;	std::exit(2);

	//testing max mem
	//int maxSiz = 1;if (verbose){		for(std::vector<char>::size_type sz = 1;   ;  sz *= 2)		{			break;			std::cerr << "attempting sz = " << sz << '\n';			std::vector<unsigned short> v(sz);		}		//cout<<"Max vec size: "<<maxSiz<<endl;	}

	if (mode == "rare_inmat"){
		Matrix* Mo = new Matrix(inF, "");//no arg for outfile &  hierachy | gene subset
		vector<DivEsts*> divvs(Mo->smplNum(),NULL);
		cout << "Using " << numThr << " ";
		vector< string > rowNames = Mo->getRowNames();

		if(rareDep == 0){
			// rarefy to smallest colSum
			rareDep = round(0.95 * Mo->getMinColSum());
			if(rareDep == 0.0){
				cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
				exit(1);
			}
		}

		// hold rarefied matrices
		// stores : repeats - sampels eg rows - vectors of columns
		int NoOfMatrices = writeFiles;
		vector< vector< map<uint, uint > > > MaRare (NoOfMatrices);
		std::vector<string> cntsNames;



		//cerr << "TH";
		std::future<rareStruct*> *tt = new std::future<rareStruct*>[numThr - 1];
		cout << "threads\n";
		uint i = 0; uint done = 0;
		while ( i < Mo->smplNum()){
			cerr << "At Sample " << i << " of " << Mo->smplNum() << " Samples";
			uint toWhere = done+numThr - 1; if ((uint)((uint)Mo->smplNum() - 2 ) < toWhere){ toWhere = Mo->smplNum() - 2; }
			for (; i < toWhere; i++){
				DivEsts * div = new DivEsts();
				tt[i - done] = async(std::launch::async, calcDivRar, i, Mo, div, rareDep, outF, repeats, writeFiles);
			}
			//use main thread to calc one sample as well
			DivEsts * div 	= new DivEsts();
			rareStruct* tmpRS;
			tmpRS 			= calcDivRar(i, Mo, div, rareDep, outF, repeats, writeFiles);


			i++;
			i = done;
			for (; i < toWhere; i++){
				rareStruct* tmpRS;
				tmpRS 		= tt[i-done].get();
				divvs[i] 		= tmpRS->div;
				string curS 	= Mo->getSampleName(i);
				divvs[i-done]->print2file(outF + curS + "_alpha_div.tsv");

				// add the matrices to the container
				if(NoOfMatrices > 0){
					if(storeBinary == true){
						binaryStoreSample(tmpMatFiles, tmpRS, rowNames,outF, cntsNames, false);
					}else{
						memoryStoreSample(tmpRS, MaRare, cntsNames, false);
					}
				}

				delete tmpRS;
			}
			// main thread divv push back
			divvs[i] = tmpRS->div;
			string curS 	= Mo->getSampleName(i);
			divvs[i]->print2file(outF + curS + "_alpha_div.tsv");

			// add the matrices to the container
			if(NoOfMatrices > 0){
				if(storeBinary == true){
					binaryStoreSample(tmpMatFiles, tmpRS, rowNames,outF, cntsNames, false);
				}else{
					memoryStoreSample(tmpRS, MaRare, cntsNames, false);
				}
			}
			delete tmpRS;
			i++;
			done = i;
		}
		printDivMat(outF + "median_alpha_diversity.tsv", divvs);
		for (size_t i = 0; i < divvs.size(); i++){
			delete divvs[i];
		}

		// write rarefaction matrices to disk
		if(NoOfMatrices > 0){
			vector< string > rowNames = Mo->getRowNames();
			if(storeBinary == true){
				printRarefactionMatrix(tmpMatFiles, outF, rareDep, cntsNames, rowNames);
			}else{
				printRarefactionMatrix(MaRare, outF, rareDep, cntsNames, rowNames);
			}
		}

		delete Mo;
		cout << "Finished\n";
		std::exit(0);
	}


	//old way of reading single samples..
	smplVec* cur = new smplVec(inF,4);
	DivEsts * div = new DivEsts();


	//placeholder for R function, not to be filled here
	std::vector<map<uint, uint>> emptyRet;
	string emptySmp;
	string skippedSample;
	cur->rarefy(rareDep,outF,repeats,div, emptyRet, emptySmp, skippedSample, writeFiles,true,false);


	div->print2file(outF+"_estimates");

	return 0;
}








void binaryStoreSample(vector< vector< string> >& tmpMatFiles, rareStruct* tmpRS, vector<string>& rowNames, string outF,  vector<string>& cntsNames, bool reshapeMap){
	// store vectors of rarefied matrix on hard disk for memory reduction
	if(reshapeMap){
		vector < string > rowIDs = tmpRS->IDs;
		vector < uint > nrowIDs(rowIDs.size());
		// convert ids into integer vector, for remapping the values
		for(uint i = 0; i < rowIDs.size(); i++){
			nrowIDs[i] = std::stoi(rowIDs[i]);
		}
		for(uint i = 0; i < tmpRS->cnts.size(); i++){
			// reshape each vector, as some are zero, and we need to rematch values and rows
			// we use the row Ids which we created correctly when splitting the vector from the input file
			std::map <uint, uint> tmpVec;
				for (auto const& x : tmpRS->cnts[i]){
					tmpVec[nrowIDs[x.first]] = x.second;
				}
				string vecLocation = printSimpleMap(tmpVec,	outF + "tmp_" + std::to_string(i) + tmpRS->cntsName + ".binary",	tmpRS->cntsName, rowNames);
				tmpMatFiles[i].push_back(vecLocation);
		}
	}else{
		for(uint i = 0; i < tmpRS->cnts.size(); i++){
			string vecLocation = printSimpleMap(tmpRS->cnts[i],	outF + "tmp_" + std::to_string(i) + tmpRS->cntsName + ".binary",	tmpRS->cntsName, rowNames);
			tmpMatFiles[i].push_back(vecLocation);
		}
	}
	// save sample name
	if(tmpRS->cntsName.size() != 0){
		cntsNames.push_back(tmpRS->cntsName);
	}
}

void memoryStoreSample(rareStruct* tmpRS, vector< vector< map<uint, uint > > >& MaRare,  vector<string>& cntsNames, bool reshapeMap){
	if(reshapeMap){
		vector < string > rowIDs = tmpRS->IDs;
		vector < uint > nrowIDs(rowIDs.size());
		// convert ids into integer vector, for remapping the values
		for(uint i = 0; i < rowIDs.size(); i++){
			nrowIDs[i] = std::stoi(rowIDs[i]);
		}
		for(uint i = 0; i < tmpRS->cnts.size(); i++){
			// reshape each vector, as some are zero, and we need to rematch values and rows
			// we use the row Ids which we created correctly when splitting the vector from the input file
			std::map <uint, uint> tmpVec;
				for (auto const& x : tmpRS->cnts[i]){
					tmpVec[nrowIDs[x.first]] = x.second;
				}
			MaRare[i].push_back(tmpVec);
		}

	}else{
		for(uint i = 0; i < tmpRS->cnts.size(); i++){
			MaRare[i].push_back(tmpRS->cnts[i]);
		}
	}
	// save sample name
	if(tmpRS->cntsName.size() != 0){
		cntsNames.push_back(tmpRS->cntsName);
	}
}



void printRarefactionMatrix(vector< vector < string > >& tmpMatFiles, string outF, int rareDep, vector<string>& cntsNames, vector<string>& rowNames){
	//vector < string > rowNames = Mo->getRowNames();
	// reassemble tmp fev files:
	for(uint i = 0; i < tmpMatFiles.size(); i++){
		string matOut = outF + "rarefied_to_" + std::to_string(rareDep) + "_n_" +  std::to_string(i) + ".tsv";
		reassembleTmpMat(tmpMatFiles[i], rowNames, cntsNames, matOut);
		// delete tmp rarefaction files now
		for(uint j = 0; j < tmpMatFiles[i].size(); j++){
			if( remove( tmpMatFiles[i][j].c_str() ) != 0 ){
				cerr << "Error deleting file: " << tmpMatFiles[i][j] << std::endl;
			}
		}
	}

}
void printRarefactionMatrix(vector<vector< map < uint, uint>>>& MaRare, string outF, int rareDep, vector<string>& cntsNames, vector<string>& rowNames){
	for(uint i = 0; i < MaRare.size(); i++){
		printRareMat(outF + "rarefied_to_" + std::to_string(rareDep) + "_n_" +  std::to_string(i) + ".tsv", MaRare[i], cntsNames, rowNames);
	}
}
