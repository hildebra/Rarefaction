#include "Rare.h"

const char* rar_ver="0.91";


rareStruct* calcDivRar(int i, Matrix* Mo, DivEsts* div, long rareDep,
	vector<vector<uint>>* abundInRow, vector<vector<uint>>* occuencesInRow,
	const vector<long>* shuffleTemplate, options* opts){
	//string outF, int repeats, int writeFiles){
	string outF ( opts->output);
	int repeats(opts->repeats);
	int writeFiles(opts->write);
	smplVec* cur = Mo->getSampleVec(i);
	string curS = Mo->getSampleName(i);
	div->SampleName = curS;
	std::vector<vector<uint>> cnts;
	vector<rare_map> cntsMap;
	string cntsName;
	string skippedNames;
	bool wrAtAll(writeFiles > 0);
	cur->rarefy(rareDep, outF, repeats,
					div, cntsMap, cntsName, skippedNames, abundInRow, occuencesInRow,
					*shuffleTemplate,
					writeFiles, false, wrAtAll );
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


rareStruct* calcDivRarVec(int i, vector<string> fileNames, DivEsts* div, long rareDep,
	vector<vector<uint>>* abundInRow, vector<vector<uint>>* occuencesInRow, 
	const vector<long>* shuffleTemplate, options* opts){
	//string outF, int repeats, int writeFiles){
	
	smplVec* cur = new smplVec(fileNames[i],4);

	std::vector<vector<uint>> cnts;
	vector< rare_map> cntsMap;
	string cntsName;
	string skippedNames;
	bool wrAtAll(writeFiles > 0);
	cur->rarefy(rareDep, outF, repeats,
					div, cntsMap, cntsName, skippedNames, abundInRow, occuencesInRow,
					writeFiles, false, wrAtAll);

	rareStruct* tmpRS 			= new rareStruct();// 	= {*div, retCnts};
	tmpRS->div 							= div;
	tmpRS->cnts 						= cntsMap;
	tmpRS->cntsName 				= cntsName;
	tmpRS->skippedNames			= skippedNames;
	tmpRS->IDs 							= cur->getRowNames();

	delete cur;

	if( remove( fileNames[i].c_str() ) != 0 ){
		cerr << "LowMem: Error deleting file: " << fileNames[i] << std::endl;
	}
	return tmpRS;
}



void helpMsglegacy(){
	string  AvailableModes = "Available run modes:\nnormalize\nsplitMat\nlineExtr\nmergeMat\nsumMat\nrarefaction\nrare_inmat\nmodule\n";
	cerr << AvailableModes << "Provide two arguments\nexiting..\n";
	cerr << "------------------------------\nAuthor: falk.hildebrand@gmail.com\n";
	std::exit(2);
}

void stateVersion(){
	printf("rarefaction tool kit (rtk) %s\n", rar_ver);
}

void helpMsg(){
	stateVersion();
	printf("\n");
	printf("USAGE\n");
	printf("    rtk <mode> -i <input.csv> -o <output> [options] \n");
	printf("\n");
	printf("MODE rarefaction\n");
	printf("\n");
	printf("OPTIONS\n");

	printf("<mode>      For rarefaction: mode can be either swap or memory.\n");
	printf("            Swap mode creates temporary files but uses less memory. \n");
	printf("            The speed of both modes is comparable.\n");
	printf("    -i      path to an .csv file to rarefy\n");
	printf("    -o      path to a output directory\n");
	printf("    -d      Depth to rarefy to. Default is 0.95 times the minimal column sum.\n");
	printf("    -r      Number of times to create diversity measures. Default is 10.\n");
	printf("    -w      Number of rarefied tables to write.\n");
	printf("    -t      Number of threads to use. Default: 1\n");
	printf("    -ns     If set, no temporary files will be used when writing rarefaction tables to disk.\n");
	printf("\n");
	printf("EXAMPLE\n");
	printf("    Rarefy a table to 1000 counts per sample with two threads. Create one table:\n");
	printf("        rtk swap -i table.csv -o outputdir/prefix. -d 1000 -r 10 -w 1 -t 2\n");
	printf("\n");
	printf("    Rarefy with most memory and least amount of IO:\n");
	printf("        rtk memory -i table.csv -o outputdir/prefix. -ns\n");
	printf("\n");
	printf("MODE: Colsums\n");
	printf("Reports the column sums of all columns in form of a sorted and an unsorted file.\n");
	printf("\n");
	printf("EXAMPLE\n");
	printf("    Repot column sums of file 'table.csv'\n");
	printf("        rtk colsums -i table.csv -o prefix\n");

	printf("\n");

	std::exit(2);
}




options::options(int argc, char** argv) :input(""), output(""), mode(""),
depth(0), repeats(10), write(0), threads(1), writeSwap(true), verbose(false),
modDB(""), modRedund(5), modEnzCompl(0.5f), modModCompl(0.5f), modWrXtraInfo(false), modCollapse(false),
xtra("") {


	bool hasErr = false;
	if (argc == 0) { return; }//could be the pseudo mod opt object

	//bool newIDthrs = false; string newIDthStr("");
	mode = argv[1];

	for (int i = 1; i < argc; i++)
	{
		if (!strcmp(argv[i], "-i"))
			input = argv[++i];
		else if (!strcmp(argv[i], "-o"))
			output = argv[++i];
		///else if (!strcmp(argv[i], "-m"))
		else if (!strcmp(argv[i], "-d"))
			depth = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-r"))
			repeats = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-w"))
			write = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-t"))
			threads = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-ns"))   // no swap
			writeSwap = false;
		else if (!strcmp(argv[i], "-v"))
			verbose = true;
		else if (!strcmp(argv[i], "-h"))
			helpMsg();
//geneMat specific args
		else if (!strcmp(argv[i], "-map"))
			map = argv[++i];
		else if (!strcmp(argv[i], "-refD"))
			referenceDir = argv[++i];
		//module specific args
		else if (!strcmp(argv[i], "-refMods"))
			modDB = (argv[++i]);
		else if (!strcmp(argv[i], "-redundancy"))
			modRedund = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-enzymeCompl"))
			modEnzCompl = (float)atof(argv[++i]);
		else if (!strcmp(argv[i], "-moduleCompl"))
			modModCompl = (float)atof(argv[++i]);
		else if (!strcmp(argv[i], "-writeExtraModEstimates"))
			modWrXtraInfo = true;
		else if (!strcmp(argv[i], "-collapseDblModules"))
			modCollapse = true;
		else if (!strcmp(argv[i], "-xtra"))
			xtra = (argv[++i]);



	}

	// sanity checks
	// we need input
	if (input == "") {//just set some defaults
		cerr << "Input must be specified\n";
		hasErr = true;
	}
	if (output == "") {//just set some defaults
		cerr << "Output must be specified\n";
		hasErr = true;
	}

	if (hasErr) {
		cerr << "Use \"rtk -h\" to get full help.\n";
		exit(98);
	}
}
void options::print_rare_details(){

	stateVersion();
	// print run mode:
	cout << "------------------------------------ "  << std::endl;
	cout << "Run information:" << std::endl;
	cout << "input file:     " << input  << std::endl;
	cout << "output file:    " << output  << std::endl;
	cout << "mode:           " << mode  << std::endl;
	if(depth != 0){
		cout << "depth:          " << depth  << std::endl;
	}else{
		cout << "depth:          0.95 x min. column sum"   << std::endl;
	}
	cout << "repeats:        " << repeats  << std::endl;
	cout << "writes:         " << write  << std::endl;
	cout << "threads:        " << threads  << std::endl;
	cout << "use swap:       ";
	if(writeSwap == false){
		cout << "false" << std::endl;
	}else{
		cout << "true" << std::endl;
	}

	//cout << "mode:           " << mode  << std::endl;
	cout << std::endl;

}



void rareExtremLowMem(string inF, string outF, int writeFiles, string arg4, int repeats, int numThr = 1, bool storeBinary = false){
	// this mode takes the file, reads it in memory
	// prints the columns to their own files
	// then it loads those files again and
	// rarefies each column
	// the measures are then combines again.

	//split mat code
	cout << "Splitting input matrix to disk" << std::endl;
	vector<string> fileNames;
	Matrix* Mo 	= new Matrix(inF, outF, "", fileNames, false, true);
	vector < string > SampleNames 	= Mo->getSampleNames();
	vector < string > rowNames 		= Mo->getRowNames();
	cout << "Done loading matrix" << std::endl;

	// abundance vectors to hold the number of occurences of genes per row
	// this will be used for ICE, ACE and or Chao2 estimation
	vector<vector<uint>> occuencesInRow(repeats, vector<uint>(Mo->rowNum(),0));
	vector<vector<uint>> abundInRow(repeats, vector<uint>(Mo->rowNum(),0));

	int rareDep 	= atoi(arg4.c_str());
	if(rareDep == 0){
		// rarefy to smallest colSum
		rareDep = (int)round(0.95f * Mo->getMinColSum());
		if(rareDep == 0.0){
			cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
			exit(1);
		}
	}
	delete Mo;


	int NoOfMatrices = writeFiles;
	vector< vector< rare_map > > MaRare (NoOfMatrices);
	std::vector<string> cntsNames;
	vector < vector < string > > tmpMatFiles (NoOfMatrices );
	int done = 0; // number of samples processed for multithreading
	uint i = 0;
	std::future<rareStruct*> *tt = new std::future<rareStruct*>[numThr - 1];


	//rarefection code
	vector<DivEsts*> divvs(fileNames.size(),NULL);
	while(i < fileNames.size()){

		// allow multithreading
		int thirds = (int) floor((fileNames.size()-3)/3);
		if(i < 3 || i % thirds == 0  ){
			cout << "At Sample " << i+1 << " of " << fileNames.size() << " Samples" << std::endl  ;
			if(i > 3 && i % thirds == 0 ){
					cout << "..." << std::endl ;
			}
		}else if( i == 3){
			cout << "..." << std::endl ;
		}

		uint toWhere = done + numThr - 1;
		if ((uint)((uint)fileNames.size() - 2 ) < toWhere){
			toWhere = fileNames.size() - 2;
		}
		// launch samples in threads
		for (; i < toWhere; i++){
			DivEsts * div 	= new DivEsts();
			div->SampleName = SampleNames[i];
			tt[i - done] = async(std::launch::async, calcDivRarVec, i, fileNames,  div, rareDep, &abundInRow, &occuencesInRow, outF, repeats, writeFiles);
		}
		// launch one in the mainthread
		DivEsts * div 	= new DivEsts();
		div->SampleName = SampleNames[i];
		rareStruct* tmpRS;
		tmpRS = calcDivRarVec(i, fileNames,  div, rareDep, &abundInRow, &occuencesInRow, outF, repeats, writeFiles);
		i++;


		// process created data, first threads, then main thread
		i = done;
		for (; i < toWhere; i++){
			rareStruct* tmpRSas;
			tmpRSas 		= tt[i-done].get();
			divvs[i] 		= tmpRSas->div;
			string curS 	= SampleNames[i];
			//divvs[i-done]->print2file(outF + curS + "_alpha_div.tsv");

			// add the matrices to the container
			if(NoOfMatrices > 0){
				if(storeBinary == true){
					binaryStoreSample(tmpMatFiles, tmpRSas, rowNames,outF, cntsNames, true);
				}else{
					memoryStoreSample(tmpRSas, MaRare, cntsNames, true);
				}
			}

			delete tmpRSas;
		}

		// main thread divv push back
		divvs[i] 			= tmpRS->div;
		string curS 	= SampleNames[i];
		//divvs[i]->print2file(outF + curS + "_alpha_div.tsv");


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

	// print the div estimates out into a file
	printDivMat(outF , divvs, true);
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

	// compute chao2, ACE, ICE and write to file
	// compute chao2, ACE, ICE and write to file
	vector<mat_fl> chao2;
	vector<mat_fl> ICE;
	vector<mat_fl> ACE;
	computeChao2(chao2, abundInRow);
	computeCE(ICE, abundInRow);
	computeCE(ACE, occuencesInRow);
	writeGlobalDiv(ICE, ACE, chao2, outF + "_gDiv.tsv");

	cout << "Finished\n";
}






int main(int argc, char* argv[])
{

	if (argc < 2) { cerr << "Not enough arguments. Use \"rtk -h\" for getting started.\n"; exit(3); }

	options* opts = new options(argc, argv);
	string inF = opts->input;
	string outF = opts->output;
	string mode = opts->mode;
	uint numThr = opts->threads;
	string arg4 = std::to_string(opts->depth);
	string map = opts->map;
	string refD = opts->referenceDir;
	//bool verbose = opts->verbose;



	//all modes that classify as rarefactions:
	if (mode == "swap" || mode == "memory") {
		opts->print_rare_details();
	}


	// start the modes
//	if (argc>3){
		if (mode == "splitMat") {
			vector<string> empt;
			Matrix* Mo = new Matrix(inF, outF, opts->xtra, empt, false);
			//Mo->splitOnHDD(outF);	//Mo->writeSums(outF);
			delete Mo;
			std::exit(0);
		//}else if (mode == "rare_lowMem") {
		//	rareLowMem(inF, outF, writeFiles,  arg4,  repeats, numThr);
	}else if (mode == "swap") {
		vector < vector < string > > tmpMatFiles(opts->write);
		rareExtremLowMem(inF, outF, opts->write,  arg4, opts->repeats, numThr, opts->writeSwap);
			std::exit(0);
		}	else if (mode == "correl2"){
			//usage: ./rare correl2 [signature matrix] [output matrix] [big gene matrix]
			//reads in signature matrix (e.g. 40 marker genes)
			//SigMatrix* Sig = new SigMatrix(inF);
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
			if (opts->modDB == "") {//try legacy mode
				Mo->estimateModuleAbund(argv, argc);// arg4, outF); //arg4 needs to be module file, outF makes the (several) matrices
			} else {
				Mo->estimateModuleAbund(opts);
			}
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
			Matrix* Mo = new Matrix(inF, outF, refD, empt, true);
			delete Mo;
			std::exit(0);
		} else if (mode == "rarefaction" || mode == "rare_inmat") {
			//rareDep = atoi(arg4.c_str());
		} else if (mode == "colSums" || mode == "colsums"  || mode == "colSum") {
		 // just load and discard the matrix and report back the colsums
		 vector<string> fileNames;
		 Matrix* Mo 	= new Matrix(inF, outF, "", fileNames, false, true, false);
		 column co 		= Mo->getMinColumn();
		 vector< pair< double, string>> colsums = Mo->getColSums();
		 Mo->writeColSums(outF);

		 cout << "" << std::endl;
		 cout << "------------------------------------" << std::endl;
		 cout << "ColSums output" << std::endl;
		 cout << "Smallest column: " << co.id << std::endl;
		 cout << "With counts:     " << co.colsum << std::endl;
		 cout << "" << std::endl;
		 cout << "Colsums where written into the files:" << std::endl;
		 cout << "    " << outF << "colSums.txt" << std::endl;
		 cout << "    " << outF << "colSums_sorted.txt" << std::endl;

		 delete Mo;
		 std::exit(0);
	 	}
		else if (mode == "geneMat"){
			cout << "Gene clustering matrix creation\n";
			if (argc < 5) {cerr << "Needs at least 4 arguments\n"; std::exit(0);}
			ClStr2Mat* cl = new ClStr2Mat(inF,outF, map, refD);
			delete cl;
			std::exit(0);
		}  else if(mode == "memory"){
		}else {
			helpMsg();
		}
//	}




	MyRNG rng;
	//test rand numbers.. check
	//std::uniform_int_distribution<unsigned long> uint_distx(0,2223951715);	for (int i=0;i<100;i++){cout<<uint_distx(rng)<<" ";}	int i;	cin>>i;	std::exit(2);

	//testing max mem
	//int maxSiz = 1;if (verbose){		for(std::vector<char>::size_type sz = 1;   ;  sz *= 2)		{			break;			std::cerr << "attempting sz = " << sz << '\n';			std::vector<unsigned short> v(sz);		}		//cout<<"Max vec size: "<<maxSiz<<endl;	}
	int rareDep = opts->depth;
	if (mode == "memory"){
		cout << "Loading input matrix to memory" << std::endl;
		Matrix* Mo = new Matrix(inF, "");//no arg for outfile &  hierachy | gene subset
		vector<DivEsts*> divvs(Mo->smplNum(),NULL);
		vector< string > rowNames = Mo->getRowNames();
		cout << "Done loading matrix" << std::endl;

		if(rareDep == 0){
			// rarefy to smallest colSum
			rareDep = (uint)round(0.95 * Mo->getMinColSum());
			if(rareDep == 0.0){
				cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
				exit(1);
			}
		}

		// hold rarefied matrices
		// stores : repeats - sampels eg rows - vectors of columns
		int NoOfMatrices = opts->write;
		vector< vector< rare_map > > MaRare (NoOfMatrices);
		std::vector<string> cntsNames;


		// abundance vectors to hold the number of occurences of genes per row
		// this will be used for Chao2 estimation
		vector<vector<uint>> abundInRow(opts->repeats, vector<uint>(Mo->rowNum(),0));
		vector<vector<uint>> occuencesInRow(opts->repeats, vector<uint>(Mo->rowNum(),0));

		//object to keep matrices
		vector < vector < string > > tmpMatFiles(opts->write);
		//cerr << "TH";
		std::future<rareStruct*> *tt = new std::future<rareStruct*>[numThr - 1];
		uint i = 0; uint done = 0;

		//get precomputed shuffle vector 
		vector<long> shuffleTemplate = Mo->suffle_pre();

		while ( i < Mo->smplNum()){
			int thirds = (int) floor(( Mo->smplNum()-3)/3);
			if(i < 3 || i % thirds == 0  ){
				cout << "At Sample " << i+1 << " of " <<  Mo->smplNum() << " Samples" << std::endl  ;
				if(i > 3 && i % thirds == 0 ){
					cout << "..." << std::endl ;
				}
			}else if( i == 3){
				cout << "..." << std::endl ;
			}
			uint toWhere = done+numThr - 1; if ((uint)((uint)Mo->smplNum() - 2 ) < toWhere){ toWhere = Mo->smplNum() - 2; }
			for (; i < toWhere; i++){
				DivEsts * div = new DivEsts();
				tt[i - done] = async(std::launch::async, calcDivRar, i, Mo, div, rareDep, &abundInRow, &occuencesInRow, &shuffleTemplate, opts);//outF, opts->repeats, opts->write);
			}
			//use main thread to calc one sample as well
			DivEsts * div 	= new DivEsts();
			rareStruct* tmpRS;
			tmpRS = calcDivRar(i, Mo, div, rareDep, &abundInRow, &occuencesInRow, &shuffleTemplate, opts); outF, opts->repeats, opts->write);


			i++;
			i = done;
			for (; i < toWhere; i++){
				rareStruct* tmpRS;
				tmpRS 		= tt[i-done].get();
				divvs[i] 		= tmpRS->div;
				string curS 	= Mo->getSampleName(i);
				//divvs[i-done]->print2file(outF + curS + "_alpha_div.tsv");

				// add the matrices to the container
				if(NoOfMatrices > 0){
					if(opts->writeSwap){
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
			//divvs[i]->print2file(outF + curS + "_alpha_div.tsv");

			// add the matrices to the container
			if(NoOfMatrices > 0){
				if(opts->writeSwap){
					binaryStoreSample(tmpMatFiles, tmpRS, rowNames,outF, cntsNames, false);
				}else{
					memoryStoreSample(tmpRS, MaRare, cntsNames, false);
				}
			}
			delete tmpRS;
			i++;
			done = i;
		}
		printDivMat(outF , divvs, true);
		for (size_t i = 0; i < divvs.size(); i++){
			delete divvs[i];
		}

		// write rarefaction matrices to disk
		if(NoOfMatrices > 0){
			vector< string > rowNames = Mo->getRowNames();
			if(opts->writeSwap){
				printRarefactionMatrix(tmpMatFiles, outF, rareDep, cntsNames, rowNames);
			}else{
				printRarefactionMatrix(MaRare, outF, rareDep, cntsNames, rowNames);
			}
		}

		delete Mo;


		// compute chao2, ACE, ICE and write to file
		vector<mat_fl> chao2;
		vector<mat_fl> ICE;
		vector<mat_fl> ACE;
		computeChao2(chao2, abundInRow);
		computeCE(ICE, abundInRow);
		computeCE(ACE, occuencesInRow);
		writeGlobalDiv(ICE, ACE, chao2, outF + "_gDiv.tsv");

		cout << "Finished\n";
		std::exit(0);
	}


	//old way of reading single samples..
	//smplVec* cur = new smplVec(inF,4);
	DivEsts * div = new DivEsts();


	//placeholder for R function, not to be filled here
	//std::vector<map<uint, uint>> emptyRet;
	//string emptySmp;
	//string skippedSample;
	//vector<vector<uint>> abundInRow;
	//cur->rarefy(rareDep,outF,repeats,div, emptyRet, emptySmp, skippedSample, &abundInRow, writeFiles,true,false);


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
			rare_map tmpVec;
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

void memoryStoreSample(rareStruct* tmpRS, vector< vector< rare_map > >& MaRare,  vector<string>& cntsNames, bool reshapeMap){
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
			rare_map tmpVec;
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
				cerr << "printRarefactionMatrix: Error deleting file: " << tmpMatFiles[i][j] << std::endl;
			}
		}
	}

}
void printRarefactionMatrix(const vector<vector< rare_map>>& MaRare, string outF, int rareDep, vector<string>& cntsNames, vector<string>& rowNames){
	for(uint i = 0; i < MaRare.size(); i++){
		printRareMat(outF + "rarefied_to_" + std::to_string(rareDep) + "_n_" +  std::to_string(i) + ".tsv", MaRare[i], cntsNames, rowNames);
	}
}
