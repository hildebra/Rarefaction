// rarefaction.cpp
// usage: rarefaction in_matrix outfile
//outfile contains richness and sample sum; in same dir several files are created
//C:\Users\Falk\SkyDrive\science\data\test\test.mat.subs

//#include "Matrix.h"
#include "Rare.h"


const char* rar_ver="0.63		 alpha";

struct cDR{
	DivEsts* div;
	std::vector<vector<vector<uint>>> retCnts;
	std::vector<string> retCntsSampleNames;
};
cDR* calcDivRar(int i, Matrix* Mo, DivEsts* div,  long rareDep, string outF,
	int repeats, int writeFiles, std::vector<vector<vector<uint>>> retCnts,
	std::vector<string> retCntsSampleNames, int NoOfMatrices){

	smplVec* cur 	= Mo->getSampleVec(i);
	string curS 	= Mo->getSampleName(i);
	div->SampleName = curS;

	cur->rarefy(rareDep, outF, repeats, div, retCnts, retCntsSampleNames,
				NoOfMatrices, writeFiles, true);

	cDR* tmpCDR 				= new cDR();// 	= {*div, retCnts};
	tmpCDR->div 				= div;
	tmpCDR->retCnts 			= retCnts;
	tmpCDR->retCntsSampleNames 	= retCntsSampleNames;

	delete cur;
	return tmpCDR;
}

/*
DivEsts* calcDivRarLegacy(int i, Matrix* Mo, DivEsts* div,  long rareDep, string outF,
	int repeats, int writeFiles, std::vector<vector<uint>> retCnts, int NoOfMatrices){
	cout << i << " ";
	smplVec* cur = Mo->getSampleVec(i);
	string curS = Mo->getSampleName(i);
	div->SampleName = curS;
	cur->rarefy(rareDep, outF, repeats, div, retCnts, NoOfMatrices, writeFiles, true);
	delete cur;
	return div;
}*/



void helpMsg(){
	string  AvailableModes = "Available run modes:\nnormalize\nsplitMat\nlineExtr\nmergeMat\nsumMat\nrarefaction\nrare_inmat\nmodule\n";
	cerr << AvailableModes << "Provide two arguments\nexiting..\n";
	cerr << "------------------------------\nAuthor: falk.hildebrand@gmail.com\n";
	std::exit(2);
}

//int main(int argc, char* argv[])
int rarefyMain(string inF, string outF, string mode,
	int repeats, long rareDep, bool verbose,
	bool returnObject, vector<vector<mat_fl>> rmatrix,
	vector< string > cnames , vector< string > rnames , DivEsts * dd,
	vector<DivEsts*> *  divvs,
	std::vector<vector<vector<uint>>> &retCnts,
	std::vector<string>& retCntsSampleNames,
	std::vector<string>& rowNames, int NoOfMatrices,
	bool transpose)
{

  /*
	//changed the interface: argv3 is now argv1 (mode),
	cout<<"Rarefaction analysis ver "<<rar_ver<<endl;
	if (argc < 3) {
		helpMsg();
	}*/

	//long rareDep = 1000;	//int repeats (20);
	//int NoOfMatrices = 1;
	/*
  //bool splitMode(false),mergeMode(false),sumUpMode(false);
	string inF = argv[2];
	string outF = argv[3];
	string mode = argv[1];
	string arg4 = "";
  */
	uint numThr = 1; //number of threads to use
	/*
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
			Mo->estimateModuleAbund(argv);// arg4, outF); //arg4 needs to be module file, outF makes the (several) matrices
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
		} else if (mode == "help" || mode == "-help" || mode == "--help"){
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
	}*/
	//int writeFiles = repeats;
	bool writeFiles = false;
	/*
  if (argc>6){
		writeFiles = atoi(argv[6]);
	}
	if (argc > 7){
		numThr = atoi(argv[7]);
	}
*/
	MyRNG rng;

  /*
  //test rand numbers.. check
	//std::uniform_int_distribution<unsigned long> uint_distx(0,2223951715);	for (int i=0;i<100;i++){cout<<uint_distx(rng)<<" ";}	int i;	cin>>i;	std::exit(2);

	//testing max mem
	//int maxSiz = 1;if (verbose){		for(std::vector<char>::size_type sz = 1;   ;  sz *= 2)		{			break;			std::cerr << "attempting sz = " << sz << '\n';			std::vector<unsigned short> v(sz);		}		//cout<<"Max vec size: "<<maxSiz<<endl;	}


*/
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
				cerr << "At Sample " << i+1 << " of " << Mo->smplNum() << " Samples\n";
			}
			uint toWhere = done+numThr - 1; if ((uint)((uint)Mo->smplNum() - 2 ) < toWhere){ toWhere = Mo->smplNum() - 2; }
			for (; i < toWhere; i++){ // with just one thread this is not used?
				DivEsts * div = new DivEsts();
				tt[i - done] = async(std::launch::async, calcDivRar, i, Mo, div, rareDep, outF,
									repeats, writeFiles, retCnts, retCntsSampleNames,
									NoOfMatrices);


			}

			//use main thread to calc one sample as well
			DivEsts * div = new DivEsts();
			//divvs[i] = calcDivRar(i, Mo, div, rareDep, outF, repeats, writeFiles);
			cDR* tmpCDr;
			//tmpCDr = new cDR;
			tmpCDr 		= calcDivRar(i, Mo, div, rareDep, outF, repeats, writeFiles, retCnts,
									retCntsSampleNames, NoOfMatrices);
			divvs->push_back(tmpCDr->div);
			retCnts 	= tmpCDr->retCnts; // Here the RAM grows and grows !is this correct or am I loosing values here?
			retCntsSampleNames = tmpCDr->retCntsSampleNames;
			delete tmpCDr;
			//cout <<'\n' <<  i << " retCnts: " << tmpCDr->retCnts.size();
			i++;
			i 			= done;
			for (; i < toWhere; i++){
				// this does not seem to be used, beause threads = 1
				(*divvs)[i] = tt[i-done].get()->div;// does this actually do anything? values are not returned to R
				string curS = Mo->getSampleName(i);
				//cout << "bob" << tt[i-done].get()->retCnts.size();
				//(*divvs)[i-done]->print2file(outF + curS + ".estimates");
			}
			i++;
			done = i;
		}
		//printDivMat(outF + "all.txt", divvs);
		/* do not delete my divvs, paul
		for (size_t i = 0; i < divvs.size(); i++){
			delete divvs[i];
		}
		*/
		//cout << "Finished\n";
		//std::exit(0);
	}




	//DivEsts * div = new DivEsts();

	// is this relevant? if this line is missing no dd is created
	//smplVec* cur = new smplVec(inF,4); // what if inF == NULL?, also does not work with matrix input
	//cur->rarefy(rareDep, outF, repeats, dd, retCnts, NoOfMatrices, writeFiles, true);
	//div->print2file(outF+"_estimates");
	return 0;
}
