#include "ClStr2Mat.h"


textBlock* getClusBlock(FILE* incl,string& lastline) {//istream& incl
	//string line;
	textBlock* ret = new textBlock;
	char buf[250];
	if (lastline.length() > 0) {
		ret->txt.push_back(lastline);
	} else {
		//getline(incl, line);
		fgets(buf, sizeof buf, incl);
		buf[strcspn(buf, "\n")] = 0;

		ret->txt.push_back(string(buf));
	}
	while (fgets(buf, sizeof buf, incl) != NULL) {
		buf[strcspn(buf, "\n")] = 0;

	//while (getline(incl, line)) {
		//if (line.substr(0, 1) == ">") {//new cluster, add line to matrix
		if(buf[0] == '>'){
			//ret->lastLine = line;
			lastline = string(buf);
			return(ret);
		}
		ret->txt.push_back(string(buf));
	}
	ret->cont = false;
	return ret;
}
void printVec(clusWrk * curClus, ofstream* mO, ofstream* gN, const vector<bool>& useSmpl) {
	long CLidx = curClus->Clnum;
	string outStr; 
	outStr.reserve(4000);
	outStr = to_string(CLidx);
	//important to have this after outStr creation!
	(*gN) << outStr + curClus->geneNamesStr;
	
	
	const vector<smat_fl>& pr = curClus->Vec;
	for (size_t i = 0; i < pr.size(); i++) {
		if (!useSmpl[i]) { continue; }
		if (pr[i] == 0) {
			outStr += "\t0";
		} else {
			outStr += "\t" + to_string(pr[i]);
		}
		//SmplSum[i] += pr[i];
	}
	outStr += "\n";
	(*mO) << outStr;
	delete curClus;
}
clusWrk* workClusBlock(textBlock* inVs, const size_t smplN,
			const string& sampleSeq,  const vector<GeneAbundance*>& GAs,
			const SmplOccurMult* smplsP, long CLidx) {
	clusWrk* ret = new clusWrk(smplN, CLidx);
	bool repFound = false;
	for (uint i = 0; i < inVs->txt.size(); i++) {
		string &line = inVs->txt[i];
		if (i == 0) {//new cluster, add headerInfo
			ret->geneNamesStr = "\t" + line;
			continue;
		}


		//1 get gene, sample (deparse)
		size_t pos = line.find("nt, >");
		size_t pos2 = line.find("...", pos + 4);
		string gene = line.substr(pos + 5, pos2 - pos - 5);

		if (!repFound && line.back() == '*') {//report representative gene
			ret->geneNamesStr += "\t" + gene + "\n";
			repFound = true;
		}

		//bool geneInAssembl(true);
		pos = gene.find(sampleSeq);
		string sample;
		//SmplOccurITmult smNum;
		pos2 = gene.find("_L", pos + 3);
		if (pos != string::npos && pos2 != string::npos) { //has the characteristic "__" sample separator
			sample = gene.substr(0, pos);
			auto smNum = smplsP->find(sample);
			if (smNum == smplsP->end()) {
#ifdef notRpackage
				cerr << "incorrect sample name: " << sample << endl << gene << endl;
				exit(55);
#endif
			}
			//2 get abundance
			//Contig + Sample Info will allow to create contig linkage between samples (CCH)
			//int contig = atoi(gene.substr(pos + 3, pos2-pos-3).c_str());
			//can be several samples (combined assembly)
			const vector<int>& smplLocs = (*smNum).second;
			for (size_t jj = 0; jj < smplLocs.size(); jj++) { //the loop takes account for multiple samples being grouped together in map
															  //CCH->addHit(smplLocs[jj], contig);
				int idxM = smplLocs[jj];
				//don't use sample, if not last in map group!
				//if (!useSmpl[idxM]) { continue; }
				smat_fl abundance = GAs[idxM]->getAbundance(gene);
				//3 add to matrix / output vector
				ret->Vec[idxM] += abundance;
				//SmplSum[idxM] += abundance;
			}
		}
		else {
			//geneInAssembl = false;
		}

	}

	
	if (!repFound) {
		ret->geneNamesStr += "\t?\n";
#ifdef notRpackage
		cerr << "No RepSeq found for " << ret->geneNamesStr << " -1" << endl;
#endif
	}
	delete inVs;
	return (ret);

}

ClStr2Mat::ClStr2Mat(options* opts):
	lastClIdWr(1),GAs(0), CCH(NULL),smplLoc(0), baseP(0), mapGr(0), SmplSum(0), smplN(0), curr(-1),
	lastline(""){
	//ifstream incl;
	FILE* incl;
	const string inF = opts->input;
	const string outF = opts->output;
	const string mapF = opts->map;
	const string basePX = opts->referenceDir;
	//set up baseP
	stringstream ss(basePX); string segments;
	while (getline(ss, segments,',') ) { baseP.push_back(segments); }


	
	matO = new ofstream(outF + ".mat");
	if (!(*matO)) {
		 #ifdef notRpackage
		cerr << "Couldn't open matrix output " << outF + ".mat" << endl;
		exit(57);
		#endif
	}
	geneNames = new ofstream(outF + ".genes2rows.txt");
	if (!(*geneNames)) {
		 #ifdef notRpackage
		cerr << "Couldn't open report file " << outF + ".genes2rows.txt" << endl;
		exit(56);
		#endif
	}
	//incl.open(inF.c_str());
	incl = fopen(inF.c_str(), "r");
	if (incl==NULL) {
		 #ifdef notRpackage
		cerr << "Couldn't open clustering file " << inF << endl;
		exit(55);
		#endif
	}
	


	//read map(s) and check that
	stringstream ss2(mapF);
	while (getline(ss2, segments, ',')) {
		read_map(segments, opts->calcCoverage, opts->calcCovMedian, opts->oldMapStyle);
	}
	this->sealMap();
	if ( baseP.size() > curr+1) {
	 #ifdef notRpackage
		cerr << "more maps than basePs\n"; 
		exit(72);
	#endif
	}


	//smplnames in out matrix
	vector<string> SmplNmsVec(smplN, "");
	for (auto it = smpls.begin(); it != smpls.end(); it++) {
		vector<int> XX = (*it).second;
		if (XX.size() == 0) { continue; }
		SmplNmsVec[ XX[XX.size()-1] ] = smplRid[(*it).first];
	}
	(*matO) << "Genes";
	for (size_t i = 0; i < SmplNmsVec.size(); i++) {
		if (!useSmpl[i]) { continue; }
		(*matO) << "\t" << SmplNmsVec[i];
	}
	(*matO) << endl;

	(*geneNames) << "#GID\tCluster\tRepSeq\n";
	CCH = new ContigCrossHit((int)smplN);
	CCH->setSmplNms(SmplNmsVec);
	//SparseMatrix * mat = new SparseMatrix();
	long CLidx = 2; //start at number 2
	const string sampleStrSep = "__";
	//vector<smat_fl> SmplSum(smplN, 0.f);
	SmplSum.resize(smplN, 0.f);
	
	int numthr = opts->threads-1;// -1 to include wrThr
	if (numthr <= 0) { numthr = 1; }

	//thread pool
	vector < job > slots(numthr);
	


	//bool repFound = false;
	textBlock* clbl;
	string lastline = "";
	//std::future <textBlock*> readThr;
	//readThr = async(std::launch::async, getClusBlock, incl, lastline);
	//readThr = async(std::launch::async, test);

	int j = 0;

	while (true) {
		clbl = getClusBlock(incl, lastline);
		//clbl = readThr.get();
		//lastline = clbl->lastLine;
		//readThr = async(std::launch::async, getClusBlock, &incl, lastline);

		
		//clusWrk* curClus = workClusBlock(clbl, smplN, CLidx,sampleStrSep,GAs, &smpls);

		//search all threads for empty slot to push job into..
		while(true) {
			if (j >= numthr) { j = 0; }
			if (slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
				slots[j].inUse = false;
				clusWrk* curClus = slots[j].fut.get();

				this->addSums(curClus);
				manage_write(curClus);
				//printVec( curClus, matO, geneNames, useSmpl);
			}
			if (slots[j].inUse == false) {
				slots[j].fut = async(std::launch::async, workClusBlock, clbl, smplN, sampleStrSep, 
					GAs, &smpls, CLidx);
				CLidx++;
				slots[j].inUse = true;
				break;//job submitted, break search for empty thread
			}
			j++;
		}
		if (!clbl->cont) { break; }
		
	}
	for (j = 0; j < numthr; j++) {//collect remaining jobs
		if (slots[j].inUse == true) {
			clusWrk* curClus = slots[j].fut.get();
			this->addSums(curClus);
			manage_write(curClus);
//			printVec(curClus, matO, geneNames, useSmpl);
			slots[j].inUse = false;
			CLidx++;
		}
	}
	//wrThr.fut.get();
	finish_write();
	

	//incl.close();
	fclose(incl);
	matO->close(); geneNames->close();
	delete matO; delete geneNames;

	ofstream* matS = new ofstream((outF + ".mat.sum"), ofstream::out);
	
	//matS.open(outF + ".mat.sum", ofstream::out);
	//print sample sum
	for (size_t i = 0; i < SmplNmsVec.size(); i++) {
		if (!useSmpl[i]) { continue; }
		//cout << SmplNmsVec[i] << "\t" << SmplSum[i] << endl;
		(*matS) << SmplNmsVec[i] << "\t" << SmplSum[i] << endl;
	}
	matS->close();  delete matS;
}
ClStr2Mat::~ClStr2Mat() {
	for (size_t i = 0; i < GAs.size(); i++) { delete GAs[i]; }
	delete CCH;
}
void ClStr2Mat::finish_write() {
	if (wrThr.inUse) {
		wrThr.fut.get();
		wrThr.inUse = false;
	}
	if (tmpSave.size() > 0) {
		auto saveCl = tmpSave.begin();
		while (saveCl != tmpSave.end()) {
			if ((*saveCl)->Clnum == (lastClIdWr + 1)) {
				printVec((*saveCl), matO, geneNames, useSmpl);
				/*wrThr.fut.get();
				wrThr.fut = async(std::launch::async, printVec, (*saveCl), matO, geneNames,
				useSmpl);
				wrThr.inUse = true;*/
				lastClIdWr++;
				tmpSave.erase(saveCl);
				saveCl = tmpSave.begin();
			}
			else {
				saveCl++;
			}
		}
	}

	if (wrThr.inUse) {
		wrThr.fut.get();
		wrThr.inUse = false;
	}
}
void ClStr2Mat::manage_write(clusWrk* curClus) {
	//push into wrThr.. (only one can run at a time, therefore wait)
	if (wrThr.inUse) {
		wrThr.fut.get();
		wrThr.inUse = false;
	}
	if (tmpSave.size() > 1) {
		auto saveCl = tmpSave.begin();
		while (saveCl != tmpSave.end()) {
			if ((*saveCl)->Clnum == (lastClIdWr + 1)) {
				printVec((*saveCl), matO, geneNames, useSmpl);
				/*wrThr.fut.get();
				wrThr.fut = async(std::launch::async, printVec, (*saveCl), matO, geneNames,
					useSmpl);
				wrThr.inUse = true;*/
				lastClIdWr++;
				tmpSave.erase(saveCl);
				saveCl = tmpSave.begin();
			}	else {
				saveCl++;
			}
		}
	}
	if (curClus->Clnum == (lastClIdWr + 1)) {
		wrThr.fut = async(std::launch::async, printVec, curClus, matO, geneNames,
			useSmpl);
		wrThr.inUse = true;
		lastClIdWr++;
	} else {
		tmpSave.push_back(curClus);
	}

}

void ClStr2Mat::read_map(const string mapF,bool calcCoverage, bool calcCovMedian, bool oldFolderStructure) {
	ifstream in;
	int map2folderIdx = 0; if (oldFolderStructure) {map2folderIdx = 1;}
	curr++;//keep track of different maps and inPaths
	uint preMapSize( (int) smplLoc.size() );

	if (curr > baseP.size()) {
		 #ifdef notRpackage
		cerr << "more maps than basePs\n";
		exit(72);
		#endif
	}
	in.open(mapF.c_str());
	if (!in) {
		 #ifdef notRpackage
		cerr << "Couldn't open mapping file " << mapF << endl;
		exit(56);
		#endif
	 }
	#ifdef notRpackage
	cout << "Reading map " << mapF << " on path " << baseP[curr] << endl;
	#endif
	SmplOccurMult CntAssGrps;
	string line; int cnt(-1); int assGrpN(-1);
	//mapping group params
	int mapGrpN(-1); //bool fillMapGrp(false); 
	SmplOccurMult CntMapGrps;
	int artiCntAssGrps(0); int skSmplCol(-1);


	while (getline(in, line)) {
		cnt ++; int sbcnt(-1);
		stringstream ss (line); string segments;
		if (line.substr(0, 1) == "#") { //read column position of smpl id, path & assGrps
			if (cnt > 0) { continue; }

			while (getline(ss, segments, '\t')) {
				sbcnt++;
				if (sbcnt==0 && segments != "#SmplID") {
					 #ifdef notRpackage
					cerr << "Map has to start with tag \"#SmplID\"\n";exit(83);
					#endif
				}
				if (sbcnt == 1 && !(segments == "Path" || segments == "SmplPrefix")) {
					#ifdef notRpackage
					cerr << "Map has to have tag \"Path\" as second entry\n";
					exit(83);
					#endif
				}
				if (segments == "AssmblGrps") {
					assGrpN = sbcnt;
                    #ifdef notRpackage
					cout << "Found Assembly groups in map\n";
                    #endif
				}
				if (segments == "MapGrps") {
					mapGrpN = sbcnt;
					//fillMapGrp = true;
                    #ifdef notRpackage
					cout << "Found Mapping groups in map\n";
                    #endif
				}
				if (segments == "ExcludeAssembly") {
					skSmplCol = sbcnt; 
					 #ifdef notRpackage
					cout << "Samples can be excluded from assembly\n";
					#endif
				}
			}
			continue;
		}


		vector<string> curLine(0);
		while (getline(ss, segments, '\t')) {
			curLine.push_back(segments);
		}
		if (skSmplCol>-1 && curLine[skSmplCol] == "1") { continue; }


		string smpID = curLine[0];


		//getline(ss, segments, '\t');
		//idx 1 for old folder structure, 0 for new folder structure
		string subDir = curLine[map2folderIdx];


		//assembly groups
		string assGrp ("");
		if (assGrpN != -1) {
			//handles assembly groups from here
			assGrp  = curLine[assGrpN];
		} else {//simulate CntAssGrps
			assGrp = itos(artiCntAssGrps);
			artiCntAssGrps++;
		}
		if (CntAssGrps.find(assGrp) != CntAssGrps.end()) {
			CntAssGrps[assGrp].push_back( (int)smplLoc.size());
		} else {
			CntAssGrps[assGrp] = vector<int>(1,(int)smplLoc.size());
		}

		if (assGrp != "" &&  CntAssGrps[assGrp].size() > 1) {
			string nsmpID = smpID + "M" + std::to_string(CntAssGrps[assGrp].size());
			if (smpls.find(nsmpID) != smpls.end()) {
#ifdef notRpackage
				cerr << "Double sample ID: " << nsmpID << endl;
				exit(12);
#endif
			}
			smpls[nsmpID] = CntAssGrps[assGrp];//(int)smplLoc.size();
			smplRid[nsmpID] = smpID;
		} else {
			if (smpls.find(smpID) != smpls.end()) {
#ifdef notRpackage
				cerr << "Double sample ID: " << smpID << endl;
				exit(12);
#endif
			}
			smpls[smpID] = vector<int>(1,(int)smplLoc.size());
			smplRid[smpID] = smpID;
		}

		//mapping groups
		string mapGrp("");
		if (mapGrpN != -1) {
			mapGrp = curLine[mapGrpN];
		}
		if (mapGrp != "" && CntMapGrps.find(mapGrp) != CntMapGrps.end()) {
			CntMapGrps[mapGrp].push_back((int)smplLoc.size());
		} else if (mapGrp != "") {
			CntMapGrps[mapGrp] = vector<int>(1, (int)smplLoc.size());
		}


		mapGr.push_back(mapGrp);
		useSmpl.push_back(true);
		smplLoc.push_back(subDir);

	}
	in.close();
	smplN = smplLoc.size();
	//read the gene abundances sample-wise in
	SmplOccur currCntMpGr;
	for (uint i = preMapSize; i < smplN; i++) {
		string pa2ab = path2counts;
		
		if (calcCoverage) { pa2ab = path2abundance; }
		if (calcCovMedian) { pa2ab = path2mediAB; }

		//only include last sample of mapping group..
		if (mapGr[i] != "") {
			SmplOccurIT cMGcnts = currCntMpGr.find(mapGr[i]);
			if (cMGcnts == currCntMpGr.end()) {
				currCntMpGr[mapGr[i]] = 1;
			}
			else {
				(*cMGcnts).second++;
			}
			if (CntMapGrps[mapGr[i]].size() != (uint) currCntMpGr[mapGr[i]]) {
				GAs.push_back(new GeneAbundance("",""));
				useSmpl[i] = false;
				continue;
			}
		}

		//read in abundance of sample
#ifdef notRpackage
		cerr << baseP[curr] + "/" + smplLoc[i] << endl;
#endif
		//DEBUG
		//continue;
		GAs.push_back(new GeneAbundance(baseP[curr] + "/" + smplLoc[i], pa2ab));
	}
}
void ClStr2Mat::sealMap() {
	SmplOccurMult nSmpls;
	for (auto smNum : smpls) {
		const vector<int>& smplLocs ( smNum.second);
		for (size_t jj = 0; jj < smplLocs.size(); jj++) { //the loop takes account for multiple samples being grouped together in map
			int idxM = smplLocs[jj];
			auto map2f = nSmpls.find(smNum.first);
			
			if (useSmpl[idxM]) {
				if (map2f == nSmpls.end()) {
					nSmpls[smNum.first] = vector<int>(1, idxM);
				}	else {
					nSmpls[smNum.first].push_back(idxM);
				}
			}
			/*else {
				nSmpls[smNum.first] = vector<int>(0);
			}*/
		}
	}
	smpls = nSmpls;
}

///////////////////////////////////////////////////
void ContigCrossHit::addHit(int Smpl, int Ctg) {

}


///////////////////////////////////////////////////

GeneAbundance::GeneAbundance(const string path, const string abunF):
	isPsAss(false){
	if (path == "" && abunF == "") { return; }//not required to read this (non-existant) file
	FILE* in;
	//first test if this is a pseudoassembly
	in = fopen((path + pseudoAssMarker).c_str(), "r");
	if (in != NULL) {
		fclose(in);
		isPsAss = true;
		return;
	}
	//not? then read abundances
	string newS = path + abunF;
	in = fopen(newS.c_str(), "r");
	if (in == NULL) {
		 #ifdef notRpackage
		cerr << "Couldn't open gene abundance file " << newS << endl;
		exit(36);
		#endif
	}
	char buf[200];
	while (fgets(buf, sizeof buf, in) != NULL) {
		//buf[strcspn(buf, "\n")] = 0;
		smat_fl abu; char gene[100];
		sscanf(buf, "%s\t%f",  gene, &abu);
		GeneAbu[gene] = abu;
		//string line(buf);
		//size_t pos = line.find("\t");
		//string gene = line.substr(0, pos);
		//GeneAbu[gene] = (smat_fl)atof(line.substr(pos + 1).c_str());
		
	}
	fclose(in);
}
smat_fl GeneAbundance::getAbundance(const string x) {
	if (isPsAss) {
		return (smat_fl) 1.f;//return one read count
	}
	SmplAbunIT fnd = GeneAbu.find(x);
	if (fnd == GeneAbu.end()) {
		return (smat_fl) 0;

 #ifdef notRpackage
cerr << "Can't find " << x << endl;
exit(33);
#endif
	}
	return (*fnd).second;
}
