// rarefaction.cpp
// usage: rarefaction in_matrix outfile
//outfile contains richness and sample sum; in same dir several files are created
//C:\Users\Falk\SkyDrive\science\data\test\test.mat.subs

//#include "Matrix.h"
#include "ClStr2Mat.h"


const char* rar_ver="0.63 alpha";


DivEsts* calcDivRar(int i, Matrix* Mo, DivEsts* div, long rareDep, string outF, int repeats, int writeFiles){
	cout << i << " ";
	smplVec* cur = Mo->getSampleVec(i);
	string curS = Mo->getSampleName(i);
	div->SampleName = curS;
	std::vector<vector<uint>> emptyRet;
	string emptySmp;
	string skippedSample;
	cur->rarefy(rareDep, outF, repeats,
					div, emptyRet, emptySmp, skippedSample,
					writeFiles, false,false);
	delete cur;
	return div;
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
	//printf("usage: rare mode input output [options] \n");
	//printf("\n");
	/*printf("Available run modes:\n");
	printf("    normalize     TODO\n");
	printf("    splitMat      TODO\n");
	printf("    lineExtr      TODO\n");
	printf("    mergeMat      TODO\n");
	printf("    sumMat        TODO\n");
	printf("    rarefaction   TODO\n");
	printf("    rare_inmat    rarefy a matrix and compute diversity measures\n");
	printf("    module        TODO\n");
	printf("\n");
	printf("Details:\n");
	printf("\n");
	printf("Mode:   rare_inmat\n");*/
	printf("usage: rare mode input output depth repeats write threads\n");
	printf("\n");
	printf("       mode      rare_inmat (Other modes are not available at the moment. Can be any string.)\n");
	printf("       input     path to a .csv file\n");
	printf("       output    path to the ouput dir and/or file prefix\n");
	printf("       depth     rarefaction depth\n");
	printf("       repeats   number of times to compute diversity\n");
	printf("       write     if the program should write the file (1 or 0)\n");
	printf("       threads   number of threads to use in parallel\n");
	printf("\n");


	std::exit(2);
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
	string mode = "rare_inmat";
	string arg4 = "";
	uint numThr = 1; //number of threads to use

	if (argc>=5){
		arg4 = argv[4];
	}

	if (argc>3){
		if (mode == "splitMat") {
			Matrix* Mo = new Matrix(inF, outF, arg4, false);
			//Mo->splitOnHDD(outF);	//Mo->writeSums(outF);
			delete Mo;
			std::exit(0);
		}
		else if (mode == "correl2"){
			//usage: ./rare correl2 [signature matrix] [output matrix] [big gene matrix]
			//reads in signature matrix (e.g. 40 marker genes)
			SigMatrix* Sig = new SigMatrix(inF);
			//Sig->estimateBinModel();
			//readMatrixLinebyLine(arg4,Sig);

		}
		else if (mode == "module") {
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
			Matrix* Mo = new Matrix(inF, outF, arg4, true);
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

	MyRNG rng;
	//test rand numbers.. check
	//std::uniform_int_distribution<unsigned long> uint_distx(0,2223951715);	for (int i=0;i<100;i++){cout<<uint_distx(rng)<<" ";}	int i;	cin>>i;	std::exit(2);

	//testing max mem
	//int maxSiz = 1;if (verbose){		for(std::vector<char>::size_type sz = 1;   ;  sz *= 2)		{			break;			std::cerr << "attempting sz = " << sz << '\n';			std::vector<unsigned short> v(sz);		}		//cout<<"Max vec size: "<<maxSiz<<endl;	}



	if (mode == "rare_inmat"){
		Matrix* Mo = new Matrix(inF, "");//no arg for outfile &  hierachy | gene subset
		vector<DivEsts*> divvs(Mo->smplNum(),NULL);
		cout << "Using " << numThr << " ";
		//cerr << "TH";
		std::future<DivEsts*> *tt = new std::future<DivEsts*>[numThr - 1];
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
			DivEsts * div = new DivEsts();
			divvs[i] = calcDivRar(i, Mo, div, rareDep, outF, repeats, writeFiles);
			i++;
			i = done;
			for (; i < toWhere; i++){
				divvs[i] = tt[i-done].get();
				string curS = Mo->getSampleName(i);
				divvs[i-done]->print2file(outF + curS + ".estimates");
			}
			i++;
			done = i;
		}
		printDivMat(outF + "all.txt", divvs);
		for (size_t i = 0; i < divvs.size(); i++){
			delete divvs[i];
		}
		cout << "Finished\n";
		std::exit(0);
	}


	//old way of reading single samples..
	smplVec* cur = new smplVec(inF,4);
	DivEsts * div = new DivEsts();


	//placeholder for R function, not to be filled here
	std::vector<vector<uint>> emptyRet;
	string emptySmp;
	string skippedSample;
	cur->rarefy(rareDep,outF,repeats,div, emptyRet, emptySmp, skippedSample, writeFiles,true,false);


	div->print2file(outF+"_estimates");

	return 0;
}
