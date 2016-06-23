#include "Matrix.h"

/*
Matrix::Matrix(const string inF):rowIDs(0),colIDs(0),sampleNameSep("")
{
	//reads matrix from HDD
	string line;
	ifstream in(inF.c_str());
	int ini_ColPerRow(0),cnt(0);

	//check MAP format
	while(getline(in,line,'\n')) {
		if(line.substr(0,1) == "#" || line.length()<2){continue;}
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss,segments,'\t')) {
			ColsPerRow++;
		}

		if (cnt==0){
			ini_ColPerRow = ColsPerRow;
		} else {
			if (ColsPerRow != ini_ColPerRow){
				cerr<<"Number of columns on line "<<cnt<<" is "<<ColsPerRow<<". Expected "<<ini_ColPerRow<<" columns.\n";
				std::exit(6);
			}
		}
		cnt++;
	}
	//vector<mat_fl> ini_vec(ini_ColPerRow-1,0.f);
	mat.resize(cnt-1,vector<mat_fl>(ini_ColPerRow-1,0.f));
	colIDs.resize(ini_ColPerRow,"");
	rowIDs.resize(cnt-1,"");
	int lineCnt= cnt;
	//reset input stream
	in.clear();
	in.seekg(0, ios::beg);
	cnt=-2;
	string segments;

	while(getline(in,line,'\n')) {
		cnt++;
		if(line.substr(0,1) == "#"){continue;}
		if (line.length()<10){continue;}
		stringstream ss;
		ss << line;
		int cnt2(-2);
		if (cnt==-1){//read header
			cnt2++;
			while (getline(ss,segments,'\t')) {
				cnt2++;
				colIDs[cnt2] = segments;
			}
			continue;
		}
		while (getline(ss,segments,'\t')) {
			cnt2++;
			if (cnt2==-1){
				rowIDs[cnt] = segments;
				continue;
			}
			mat[cnt][cnt2] = (mat_fl)atof(segments.c_str());
		}

	}
	in.close();

}
*/

inline mat_fl median(std::vector<mat_fl> vec, bool ignoreZeros)
{
	if (vec.size() == 0) { return (mat_fl) 0; }
	//std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, vec.end());
	sort(vec.begin(), vec.end());
	size_t i (0);
	if (ignoreZeros) {
		for (; i < vec.size(); i++) {
			if (vec[i] > 0) {
				break;
			}
		}
		if (vec.size() == i) { return (mat_fl)0; }
	}
	size_t size = vec.size() - i;

	if (size % 2 == 0) {
		return (vec[(size / 2) - 1 + i] + vec[size / 2 + i]) / 2;
	}
	return  vec[size / 2 + i];
}
void vecPurge(vector<vector<mat_fl>>& vec, mat_fl val) {
	for (size_t i = 0; i < vec.size(); i++) {
		for (size_t j = 0; j < vec[i].size(); j++) {
			vec[i][j] -= val;
		}
	}
}
string join(const vector<string>& in, const string &X) {
	string ret(in[0]); for (size_t i = 1; i < in.size(); i++) { ret += X + in[i]; } return ret;
}



//*********************************************************
ModStep::ModStep(const string & s):alternates(0), redundancy(0) {
	istringstream ss(s);
	string token(""),tok2("");
	//std::istream_iterator<std::string> beg(ss), end;
	//std::vector<std::string> tokens(beg, end); // done!
	//for (auto& s : tokens) { cout << s; }
	while (std::getline(ss, token, '\t')) {
		stringstream buff(token);
		vector<string> tmp(0);
		while (getline(buff, tok2, ',')) {
			tmp.push_back(tok2);
		}
		alternates.push_back(tmp);
	}
}

//how often does the KO used occur in total dataset?
void ModStep::setRedund(ModOccur& m) {
	redundancy.resize(alternates.size());
	for (size_t i = 0; i < alternates.size(); i++) {
		vector<int> tmp(alternates[i].size(), 0);
		for (size_t j = 0; j < alternates[i].size(); j++) {
			tmp[j] = m[alternates[i][j]];
		}
		redundancy[i] = tmp;
	}
}
void ModStep::getAllKOs(list<string>& ret) {

	for (size_t i = 0; i < alternates.size(); i++) {
		for (size_t j = 0; j < alternates[i].size(); j++) {
			ret.push_back(alternates[i][j]);
		}
	}
}
void ModStep::abundParts(const vector<mat_fl>& v, const unordered_map<string, int>& IDX,
	vector<mat_fl>& abund, vector<bool>& active, vector<string>& KOdescr,
	float hitComplRatio, int redund) {
	//some params, should be fine tuned if possible
	//float hitComplRatio(0.8f);

	active.resize(alternates.size(), false);
	abund.resize(alternates.size(), (mat_fl)0);
	//actual deep routine to determine if KOs in this step satisfy presence conditions
	for (size_t i = 0; i < alternates.size(); i++) {
		size_t altS = alternates[i].size(); float hits(0);
		vector<mat_fl> tmpAB(altS, (mat_fl)0);
		for (size_t j = 0; j < altS; j++) {
			//redundant???
			if (redundancy[i][j] > redund) {
				altS--; //reduce size of this set

				continue;
			}
			//check for actual abundance in matrix-vector subpart
			auto fn = IDX.find(alternates[i][j]);
			if (fn == IDX.end()) {
				tmpAB[j] = 0;
			} else {
				hits++;
				tmpAB[j] = v[fn->second];
			}
		}
		if (altS == 0) {
			//r[i] = (mat_fl)-1;//signal that removed due to redundancy
			continue;
		}
/*		if (hits > 0) {
			int x = 0;
		}*/
		if ( hits / (float)altS >= hitComplRatio) {
			abund[i] = median(tmpAB);
			active[i] = true;
		}
	}

	//return (active);
}

//*********************************************************
Module::Module(vector<string>& n) :name(""), description(""), steps(0){
	string token("");
	for (size_t i = 0; i < n.size(); i++) {
		if (i == 0) {//module name & description
			istringstream ss(n[i]);
			getline(ss, token, '\t'); name = token;
			getline(ss, token, '\t'); description = token;
		} else { //actual modules components (e.g. KOs)
			steps.push_back( ModStep(n[i] ));
		}
	}
}
mat_fl Module::pathAbundance(const vector<mat_fl>& v,  const unordered_map<string, int>& IDX,
		const int redund, const float PathwCompl, const float enzymCompl, string & savedNmsKO, float& modScoreOut) {

	//initial parameters
	//float PathwCompl(0.6f); //corresponds to -c
	//float enzymCompl(0.8f);
	//int redund(0);

	vector< vector< mat_fl >> abunds(steps.size(), vector<mat_fl>(0));//contains abundance
	vector< vector< bool >> active(steps.size());//contains info if path was even active
	vector<mat_fl> preMed(steps.size(), (mat_fl)0), postMed(steps.size(), (mat_fl)0);
	//auto t = IDX.find("xx");
	for (size_t i = 0; i < steps.size(); i++) {
		steps[i].abundParts(v, IDX,abunds[i], active[i], altKOs[i], enzymCompl,  redund);
		//determine median overall value
		preMed[i] = median(abunds[i]);
	}
	mat_fl pm = median(preMed);
	mat_fl retval(0);


	if (0) {
		//VAR 1
		//select one median value per step only, for the final pathway median abundance
		for (size_t i = 0; i < steps.size(); i++) {
			vector<mat_fl> tmp(abunds[i].size(), (mat_fl)0);
			for (size_t j = 0; j < abunds[i].size(); j++) {
				tmp[j] = abunds[i][j] - pm;
				exit(99);
			}
		}
	}
	else {
		//VAR 2
		//calc abundance for each possible (active) path and substract from others, but add remainders up
		//vector<mat_fl> tmp(steps.size(), (mat_fl)0);
		vector<mat_fl> curP(steps.size(),(mat_fl)0);
		float act(0.f),shldAct((float)steps.size());
		vector<int> decIdx(steps.size(), 0);
//ini decIdx
		for (size_t i = 0; i < steps.size(); i++) {
			size_t dI = 0; double maxAB = 0;
			while (dI < int (active[i].size() )) {
				if (abunds[i][dI] > maxAB && active[i][dI]) {
					decIdx[i] = dI;
					maxAB = abunds[i][decIdx[i]];
				}
				dI++;
			}
		}
		bool saveKOnames(true);// save names of KOs used in extra file?

		while (1) { // this loop goes over every possible path combination
			for (size_t i = 0; i < steps.size(); i++) {
				if (abunds[i][decIdx[i]] > 0 && active[i][decIdx[i]]) {
					curP [i] = abunds[i][ decIdx[i] ]; act++;
					if (saveKOnames) {
						savedNmsKO += altKOs[i][decIdx[i]];// +",";
					}
				} /*else if (!active[i][decIdx[i]]) {
					shldAct--;
				} */
			}
			if ( (act / shldAct) >= PathwCompl) { //shldAct > 0 &&
				//this part parses out the KOs that are actually active
				mat_fl curM = median(curP,true);
				retval += curM;
				vecPurge(abunds,curM);
				modScoreOut = act / shldAct;
			}
			else {
				savedNmsKO = "";
				modScoreOut = 0.f;
			}
			break;
		}

	}
	return retval;
}


//*********************************************************
//read in module file
Modules::Modules(const string& inF):
	redund(1),PathwCompl(0.6f), enzymCompl(0.8f) {
	ifstream is(inF.c_str());
	string line(""); vector<string> buffer(0);
	string ModToken = "M";
	while (safeGetline(is, line)) {
		//comment
		if (line[0] == '#') { continue; }
		//new module opens, create old module
		if (line.find(ModToken) == 0) {
			if (buffer.size() > 0) {
				mods.push_back(Module(buffer));
			}
			buffer.resize(0);
		}
		if (line.size() > 3) {
			buffer.push_back(line);
		}
	}
	//create last module
	mods.push_back(Module(buffer));
	buffer.resize(0);
	is.close();
	//redundancy of KOs
	calc_redund();

	//set up names
	moduleNames.resize(mods.size(), "");
	for (size_t i = 0; i < mods.size(); i++) {
		moduleNames[i] = mods[i].name;
	}
	//set up descriptions
	moduleDescriptions.resize(mods.size(), "");
	for (size_t i = 0; i < mods.size(); i++) {
		moduleDescriptions[i] = mods[i].description;
	}
}
void Modules::calc_redund() {
	list<string> fL; //full list of KOs
	for (size_t i = 0; i < mods.size(); i++) {
		mods[i].getAllKOs(fL);
	}
	MO.clear();
	//calc abundance
	for (auto& s : fL) {
		auto fnd = MO.find(s);
		if (fnd == MO.end()) {
			MO[s] = 1;
		}
		else {
			fnd->second++;
		}
	}
	//stats on KO redundancy (print out only)
	vector<int>statKOr = vector<int>(0, 0);
	int maxRed = 0;
	for (auto kor : MO) {
		if (kor.second > maxRed) {
			statKOr.resize((kor.second + 1), 0);
			maxRed = kor.second;
		}
		statKOr[kor.second] ++;
	}
	//cout << "stats on DB KO redundancy (redundancy : occurence):\n";
	for (size_t i = 0; i < statKOr.size(); i++) {
		if (statKOr[i] < 1) { continue; }
		//cout << i << " : " << statKOr[i] << endl;
	}
	for (size_t i = 0; i < mods.size(); i++) {
		mods[i].setReddundancy(MO);
	}

}

vector<mat_fl> Modules::calcModAbund(const vector<mat_fl>& v, const unordered_map<string,
	int>& IDX, vector<string> &retStr, vector<float> &retScore) {
	vector<mat_fl> ret(mods.size(),(mat_fl)0);
	retStr.resize(mods.size(), ""); retScore.resize(mods.size(), 0.f);
	mat_fl unass_cnt(0.f);//TOGO
	//vector<bool> usedKOs(v.size()); //initial idea to save KOs used to estimated unassigned fraction - better to scale by seq depth external

	for (size_t i = 0; i < mods.size(); i++) {
		ret[i] = mods[i].pathAbundance(v,IDX, redund, PathwCompl, enzymCompl, retStr[i], retScore[i]);
	}
	return ret;
}


//*********************************************************
Matrix::Matrix(void)
	:rowIDs(0), colIDs(0), maxCols(0), HI(0), maxLvl(0), sampleNameSep(""), doSubsets(false), doHigh(false)
{
}

Matrix::Matrix(const string inF, const string outF, const string xtra, bool highLvl)
	: rowIDs(0), colIDs(0), maxCols(0), HI(0), maxLvl(0), sampleNameSep(""), doSubsets(false), doHigh(highLvl)
{
	//reads matrix from HDD
	//and writes it simultaneously to single files
	if (doHigh){
		read_hierachy(xtra);
	} else if (xtra.length() > 2){
		read_subset_genes(xtra);
	}
	string line;
	ifstream in(inF.c_str());
	if (!in){ cerr << "Cant open file " << inF << endl; std::exit(11); }
	int ini_ColPerRow(0),cnt(0);


	//check MAP format
	//while (safeGetline(in, line)) {
	while (getline(in, line,'\n')) {
		if (line.substr(0, 1) == "#" || line.length()<2){ continue; }
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss,segments,'\t')) {
			ColsPerRow++;
		}

		if (cnt==0){
			ini_ColPerRow = ColsPerRow;
		} else {
			if (ColsPerRow != ini_ColPerRow){
				cerr<<"C1: Number of columns on line "<<cnt<<" is "<<ColsPerRow<<". Expected "<<ini_ColPerRow<<" columns.\n"<<line<<endl;
				std::exit(6);
			}
		}
		cnt++;
		if (cnt>10){break;}
	}
	colIDs.resize(ini_ColPerRow-1,"");
	vector<double> colSum(ini_ColPerRow-1,0.0);
	vector<ofstream> outFs(ini_ColPerRow-1);
	vector<string> outStr(ini_ColPerRow-1);
	//int lineCnt= cnt;
	//reset input stream
	in.clear();
	in.seekg(0, ios::beg);
	cnt=-1;
	string segments;
	safeGetline(in, line);
	//getline(in, line, '\n');
	while (line.substr(0, 1) == "#"){
		safeGetline(in, line);
	}
	stringstream sso;
	int cnt2(-2);
	sso << line;
	while (getline(sso,segments,'\t')) {
		cnt2++;
		if (segments.length() > 150){
			cerr << segments << " error!\n"; std::exit(5);
		}
		if (cnt2==-1){continue;}
		colIDs[cnt2] = segments;
		string oF2 = outF + sampleNameSep + colIDs[cnt2];
		if (!doHigh){
			outFs[cnt2].open(oF2.c_str(), ios_base::out);
			outFs[cnt2].precision(12);
		}
	}
	if (doHigh){
		for (int i = 0; i < maxLvl; i++){
			HI.push_back(new HMat(LvlNms[i], colIDs, vector<string>(0)));
		}
	}
	string rowID="";
	int geneCnt(0);
	int cntNA(0);
	while (safeGetline(in, line)) {
	//while (getline(in, line, '\n')) {
		cnt++;
		if(line.substr(0,1) == "#"){continue;}
		if (line.length()<10){continue;}
		int cnt2(-2);
		vector<string> taxa(0);
		stringstream ss;
		ss << line;
		bool breaker(false);
		std::map<std::string, vector<string>>::iterator fnd;
		while (getline(ss, segments, '\t')) {
			cnt2++;
			if (cnt2 == -1){
				rowID = segments;
				if (doSubsets && subset.find(rowID) == subset.end()){
					breaker = true;
					break;
				}
				if (doHigh){
					fnd = LUp.find(rowID);
					if (fnd == LUp.end()){//needs to be added to HMat
						taxa = vector<string>(maxLvl, "-1");
						cntNA++;
						if (cntNA < 100) {
							cout << "Row ID " << rowID << " is not in hierachy.\n";// \nAborting..\n"; std::exit(24);
							if (cntNA == 99) { cout << " ..\n"; }
						}
					}
					else {
						taxa = (*fnd).second;
					}
					/*
					string lngTax = "";
					for (int tt = 0; tt < maxLvl; tt++){
						lngTax += taxa[tt] ;
						taxa[tt] = lngTax;
						lngTax += +";";
					}
					*/
				}
				geneCnt++;
				continue;
			}
			mat_fl tmp =  atof(segments.c_str());
			if (doHigh){//1:finds relevant rowID, extracts taxa; 2:add on all HighLvl mats
				for (int tt = 0; tt< maxLvl; tt++){
					HI[tt]->set(taxa[tt], cnt2, tmp);
				}
				colSum[cnt2] += (double)tmp;
			}
			else if (tmp>0){//write to File
				//outFs[cnt2]<<rowID<<"\t"<<tmp<<endl;
				outStr[cnt2] += rowID+"\t"+segments.c_str()+"\n";
				colSum[cnt2] += (double)tmp;
			}
		}
		if (breaker){
			continue;
		}
		if (cnt2+2 != ini_ColPerRow){
			cerr<<"C2: Number of columns on line "<<cnt<<" is "<<cnt2+2<<". Expected "<<ini_ColPerRow<<" columns.\n";
			std::exit(62);
		}
		if (cnt % 10000 ==0){
			for (size_t cnt2=0;cnt2<outStr.size();cnt2++){
				outFs[cnt2] << outStr[cnt2];outStr[cnt2] = "";
			}
		}

	}
	in.close();
	ofstream out;
	if (doHigh){//write out high lvl mats
		for (int i = 0; i < maxLvl; i++){
			string oF2 = outF + LvlNms[i] + ".txt";
			out.open(oF2.c_str(), ios_base::out);
			HI[i]->print(out);
			out.close();
		}
	}
	else {//close filestreams to single sample files
		for (size_t i = 0; i < outFs.size(); i++){
			outFs[i] << outStr[i];
			outFs[i].close();
		}
	}
	//write colSums
	string oF2 = outF + sampleNameSep + "sums.txt";

	out.open(oF2.c_str(),ios_base::out);
	out.precision(12);
	for (size_t smpl=0;smpl<(colIDs.size()); smpl++){
		out<<colIDs[smpl]<<"\t"<<colSum[smpl]<<endl;
	}
	out.close();
	cout << "Read " << geneCnt << " genes" << endl;
}

Matrix::Matrix(const vector<string>& rnms, const vector<string>& cnms):
	rowIDs(rnms),colIDs(cnms), maxCols((int)cnms.size())
{
	vector<mat_fl> iniV = vector<mat_fl>(rnms.size(), (mat_fl)0);
	mat.resize(maxCols, iniV);
}
Matrix::Matrix(const string inF, const string xtra, bool highLvl)
	: rowIDs(0), colIDs(0), maxCols(0), HI(0), maxLvl(0), sampleNameSep(""), doSubsets(false), doHigh(highLvl)
{
	//reads matrix from HDD
	//and writes it simultaneously to single files
	if (doHigh){
		read_hierachy(xtra);
	}
	else if (xtra.length() > 2){
		read_subset_genes(xtra);
	}
	string line;
	ifstream in(inF.c_str());
	if (!in){ cerr << "Cant open file " << inF << endl; std::exit(11); }
	int ini_ColPerRow(0), cnt(0);


	//check MAP format
	//while (safeGetline(in, line)) {
	while (getline(in, line, '\n')) {
		if (line.substr(0, 1) == "#" || line.length()<2){ continue; }
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss, segments, '\t')) {
			ColsPerRow++;
		}

		if (cnt == 0){
			ini_ColPerRow = ColsPerRow;
		}
		else {
			if (ColsPerRow != ini_ColPerRow){
				cerr << "C1: Number of columns on line " << cnt << " is " << ColsPerRow << ". Expected " << ini_ColPerRow << " columns.\n" << line << endl;
				std::exit(63);
			}
		}
		cnt++;
		if (cnt>10){ break; }
	}
	if (ini_ColPerRow == 0) {
		cerr << "Empty matrix provided\n";
		return;
	}
	colIDs.resize(ini_ColPerRow - 1, "");
	vector<double> colSum(ini_ColPerRow - 1, 0.0);
	in.clear();
	in.seekg(0, ios::beg);
	cnt = -1;
	string segments;
	safeGetline(in, line);
	//getline(in, line, '\n');
	while (line.substr(0, 1) == "#"){
		safeGetline(in, line);
	}
	stringstream sso;
	int cnt2(-2);
	//read & prep header
	sso << line;
	while (getline(sso, segments, '\t')) {
		cnt2++;
		if (segments.length() > 150){
			cerr << segments << " error!\n"; std::exit(5);
		}
		if (cnt2 == -1){ continue; }
		colIDs[cnt2] = segments;
	}
	if (doHigh){
		for (int i = 0; i < maxLvl; i++){
			HI.push_back(new HMat(LvlNms[i], colIDs, vector<string>(0)));
		}
	}
	string rowID = "";
	int geneCnt(0);
	int cntNA(0);
	//vector<mat_fl> emptyVec(ini_ColPerRow, (mat_fl)0);
	mat.resize(ini_ColPerRow -1, vector<mat_fl>(0));
	while (safeGetline(in, line)) {
		//while (getline(in, line, '\n')) {
		cnt++;
		if (line.substr(0, 1) == "#"){ continue; }
		if (line.length()<5){ continue; }
		int cnt2(-2);
		vector<string> taxa(0);
		stringstream ss;
		ss << line;
		bool breaker(false);
		std::map<std::string, vector<string>>::iterator fnd;
		while (getline(ss, segments, '\t')) {
			cnt2++;
			if (cnt2 == -1){
				rowID = segments;
				if (rowID == "mapped"){//to deal with mocat files
					breaker = true;
					break;
				}
				if (doSubsets && subset.find(rowID) == subset.end()){
					breaker = true;
					break;
				}
				//set index for rows..
				rowID_hash[rowID] = cnt;
				rowIDs.push_back(rowID);

				if (doHigh){
					fnd = LUp.find(rowID);
					if (fnd == LUp.end()){//needs to be added to HMat
						taxa = vector<string>(maxLvl, "-1");
						cntNA++;
						if (cntNA < 100) {
							cout << "Row ID " << rowID << " is not in hierachy.\n";// \nAborting..\n"; std::exit(24);
							if (cntNA == 99) { cout << " ..\n"; }
						}
					}
					else {
						taxa = (*fnd).second;
					}
					/*
					string lngTax = "";
					for (int tt = 0; tt < maxLvl; tt++){
					lngTax += taxa[tt] ;
					taxa[tt] = lngTax;
					lngTax += +";";
					}
					*/
				}
				geneCnt++;
				continue;
			}
			mat_fl tmp = atof(segments.c_str());
			if (doHigh){//1:finds relevant rowID, extracts taxa; 2:add on all HighLvl mats
				for (int tt = 0; tt< maxLvl; tt++){
					HI[tt]->set(taxa[tt], cnt2, tmp);
				}
				colSum[cnt2] += (double)tmp;
			}
			else {//would make it a sparse matrix: if (tmp>0) (but looses indexing)
				mat[cnt2].push_back(tmp);
				colSum[cnt2] += (double)tmp;
			}
		}

		if (breaker){
			continue;
		}
		if (cnt2 + 2 != ini_ColPerRow){
			cerr << "C2: Number of columns on line " << cnt << " is " << cnt2 + 2 << ". Expected " << ini_ColPerRow << " columns.\n";
			std::exit(64);
		}

	}
	in.close();
	maxCols = (int)mat.size();
	cout << "Read " << geneCnt << " genes" << endl;

	if (geneCnt == 0) {
		cerr << "No genes read.. aborting\n";
		exit(0);
	}
}

void Matrix::estimateModuleAbund(char ** argv, int argc) {
	string moduleFile = argv[4];
	string outFile = argv[3];
	string doModKOest = argv[4];
	//read module DB
	Modules* modDB = new Modules(moduleFile);

	//ini options
	int redundancy = atoi(argv[5]);
	float pathCompl = (float)atof(argv[6]);
	float enzyCompl = (float)atof(argv[7]);
	bool writeKOused = false;
	//cout << argv[8] << endl;
	if (argc > 8 && strcmp(argv[8],"1")==0 ) { writeKOused = true; }

	modDB->setEnzymCompl(enzyCompl);
	modDB->setPathwCompl(pathCompl);
	modDB->setRedund(redundancy);
	//matrix vector of module abudnance
	//vector<vector<mat_fl>> modMat(maxCols);
	Matrix modMat = Matrix(modDB->modNms(), colIDs);
	vector<vector<string>> modStr(maxCols);
	vector<vector<float>> modScore(maxCols);
	for (int i = 0; i < maxCols; i++) {
//		cerr << i << " ";
		//TODO: add unknown counts
		modMat.addTtlSmpl(
			modDB->calcModAbund(mat[i], rowID_hash, modStr[i], modScore[i])

			,i );
	}
	//write description
	ofstream of; ofstream of2; vector<string> moD = modDB->modDescr(); vector<string> moN = modDB->modNms();
	string nos = outFile+".descr";
	of.open(nos.c_str());
	for (size_t i = 0; i < moD.size(); i++) {
		of << moN[i] << "\t" << moD[i] << endl;
	}
	of.close();
	//write KOs used
	if (writeKOused) {
		nos = outFile + ".KOused";
		of.open(nos.c_str());
		nos = outFile + ".MODscore";
		of2.open(nos.c_str());

		//write SmplIDs
		for (size_t i = 0; i < colIDs.size(); i++) {
			of << "\t" << colIDs[i]; of2 << "\t" << colIDs[i];
		}
		of << endl; of2 << endl;

		for (size_t i = 0; i < moD.size(); i++) {
			bool hasKOUse = false;
			for (size_t k = 0; k < modStr.size(); k++) {
				if (modStr[k][i] != "") { hasKOUse=true; break; }
			}
			if (!hasKOUse) { continue; }
			of << moN[i]; of2 << moN[i];
			for (size_t k = 0; k < modStr.size(); k++) {
				of << "\t" << modStr[k][i] ;
				of2 << "\t" << modScore[k][i] ;
			}
			of << endl; of2 << endl;
		}
		of.close(); of2.close();
	}

	//write module matrix
	cerr << "Write Matrix\n";
	modMat.writeMatrix(outFile+".mat",true);
	delete modDB;

}
void Matrix::addColumn(string cname) {
	//increase maxLvl cnt
	maxCols++;
	colIDs.push_back(cname);
	colID_hash[cname] = maxCols-1;
	for (uint i = 0; i < mat.size(); i++) {
		//mat[i][smpl];
		mat[i].resize(maxCols, 0);
	}
}

Matrix::~Matrix(void)
{
	for (unsigned int i = 0; i < HI.size(); i++){
		delete HI[i];
	}

}
void Matrix::addRow(vector<mat_fl> x) {
	mat.push_back(x);
	//maybe add some security checks later..
}

void Matrix::normalize() {
	vector<mat_fl> allSums(colIDs.size(), (mat_fl)0);
	for (size_t smpl = 0; smpl<(colIDs.size() - 1); smpl++) {
		mat_fl sums(0);
		for (size_t i = 0; i<rowIDs.size(); i++) {
			sums += mat[smpl][i];
		}
		allSums[smpl] = sums;
	}
	for (size_t smpl = 0; smpl < (colIDs.size() - 1); smpl++) {
		for (size_t i = 0; i < rowIDs.size(); i++) {
			mat[smpl][i] /= allSums[smpl];
		}
	}
	/* //DEBUG
	for (size_t smpl = 0; smpl < (colIDs.size() - 1); smpl++) {
		mat_fl sums(0);
		for (size_t i = 0; i < rowIDs.size(); i++) {
			sums += mat[smpl][i];
		}
		sums++;
	}
	*/
}

void Matrix::transpose(){
	// takes the matrix and transposes it
	// column ID and row ID have to be swapped as well
	vector< vector< mat_fl > >  transpMat(mat[0].size(), vector< mat_fl >(mat.size()));

	for(int i = 0; i < mat.size(); i++){
		for(int j = 0; j < mat[i].size(); j++){
			transpMat[j][i] = mat[i][j];
		}
	}

	// switch column and row names
	vector< string >  rowIDst 	= rowIDs;
	rowIDs 						= colIDs;
	colIDs						= rowIDst;

	// swap the matrices
	mat = transpMat;
}

vector<mat_fl> Matrix::getRowSums() {
	vector<mat_fl> rowSums(rowIDs.size(), 0);
	for (size_t i = 0; i < rowIDs.size(); i++) {
		for (size_t smpl = 0; smpl < (colIDs.size() ); smpl++) {
			rowSums[i] += mat[smpl][i];
		}
	}
	return rowSums;
}

void Matrix::writeMatrix(const string of, bool onlyFilled) {
	ofstream out;
	out.open(of.c_str(), ios_base::out);
	out.precision(8); out << "Gene";
	for (size_t smpl = 0; smpl < (colIDs.size() ); smpl++) {
		out << "\t" << colIDs[smpl ];
	}
	out << endl;
	vector<mat_fl> rowSums;
	size_t cidS(colIDs.size());
	if (onlyFilled) { rowSums = getRowSums(); }
	for (size_t i = 0; i<rowIDs.size(); i++) {
		if (onlyFilled && rowSums[i]==0) {
			continue;
		}
		out << rowIDs[i] ;
		for (size_t smpl = 0; smpl<cidS; smpl++) {
			out << "\t" << mat[smpl][i];
		}
		out << endl;
	}
	out.close();
}

void Matrix::read_subset_genes(const string xtra){
	string line;
	ifstream in(xtra.c_str());
	string ID = ""; int cnt = 0;
	if (!in){ cerr << "Can't open geneID file " << xtra << endl; std::exit(13); }
	while (getline(in, ID, '\n')) {
		subset[ID] = 1;
		cnt++;
	}
	in.close();
	//just check that something was read
	if (cnt >= 1){
		doSubsets = true;
		cout << "Read " << cnt << " gene subsets\n";
	}
}
void Matrix::read_hierachy(const string xtra){
	int maxHir = 7;
	vector<string> features;
	string line;
	ifstream in(xtra.c_str());
	int cnt = 0;
	if (!in){ cerr << "Can't open hierachy file " << xtra << endl; std::exit(13); }
	for (int k = 0; k < maxHir; k++){
		//string tmp = ;
		LvlNms.push_back("L" + stringify((double)k));
	}

	while (getline(in, line, '\n')) {
		if (line.substr(0, 1) == "#"){ continue; }
		vector<string> pseudo(maxHir, "?");
		cnt++;
		string segs;
		string segs2;
		stringstream ss;
		ss << line;
		getline(ss, segs, '\t');
		getline(ss, segs2, '\t');
		string spl;
		stringstream hir; hir << segs2;
		int i = 0;
		//string lngTax="";
		while (getline(hir, spl, ';')){
			//add previous levels
			pseudo[i] = spl;
			i++;
			//lngTax += spl + ";";
			//index possible features
			if (i >= maxHir){ break; }
		}
		if (i > maxLvl){ maxLvl = i; }
		LUp[segs] = pseudo;
	}
	in.close();
	cout << "Read hierachy. Found " << maxLvl << " hierachical levels.\n";
}
void Matrix::splitOnHDD(string out_seed){
	for (size_t smpl=0;smpl<(colIDs.size()-1); smpl++){
		string oF2 = out_seed + sampleNameSep + colIDs[smpl+1];
		ofstream out;
		out.open(oF2.c_str(),ios_base::out);
		out.precision(12);
		for (size_t i=0;i<rowIDs.size();i++){
			if (mat[i][smpl]==0){continue;}
			out<< rowIDs[i]<<"\t"<<mat[i][smpl]<<endl;
		}
		out.close();
		//if (smpl>20){std::exit(9);}
	}
}


void Matrix::writeSums(string out_seed){
	string oF2 = out_seed + sampleNameSep + "sums.txt";
	ofstream out;
	out.open(oF2.c_str(),ios_base::out);
	out.precision(12);
	//for (OTUid = colID_hash.begin(); OTUid != colID_hash.end(); OTUid++) {
	for (size_t smpl = 0; smpl<(colIDs.size() - 1); smpl++) {
		mat_fl sums(0);
		for (size_t i=0;i<rowIDs.size();i++){
			sums+=mat[i][smpl];
		}
		out<<colIDs[smpl+1]<<"\t"<<sums<<endl;
	}
	out.close();
}


//Saprse Matrix class
SparseMatrix::SparseMatrix() :mat(0), colNames(0), rowIDs(0){}

void SparseMatrix::addCount(string smpl, int row, smat_fl abund) {
	//find sample
	SmplAbunIT tar = mat[row].find(smpl);
	if (tar == mat[row].end()) {//create entry & expand matrix
		mat[row][smpl] = abund;
	} else {
		(*tar).second += abund;
	}

//keep track of listed samples
	SmplOccurIT tar2 = colNames.find(smpl);
	if (tar2 == colNames.end()) {
		colNames[smpl] = 1;
	} else {
		(*tar2).second++;
	}
}


HMat::HMat(string L, vector<string> Samples, vector<string> Features)
:LvlName(L), FeatureNs(Features), SampleNs(Samples),mat(0){
	empty = vector<mat_fl>(SampleNs.size(), 0);
	mat.resize(FeatureNs.size(), empty);
	for (unsigned int i = 0; i < FeatureNs.size(); i++){
		Feat2mat[FeatureNs[i]] = i;
	}
}

void HMat::set(string kk, int j, mat_fl v) {
	mat_fl div(1); size_t pos(kk.find(",", 0)), npos(0);
	vector<string> subkk(0);
	while (pos != string::npos ){
		subkk.push_back(kk.substr(npos, pos-npos));
		npos = pos+1;
		pos = kk.find(",", npos);
		div += 1.f;
	}
	subkk.push_back(kk.substr(npos));

	for (uint t = 0; t < subkk.size(); t++){
		string yy = subkk[t];
		std::map<string, int>::iterator i = Feat2mat.find(yy);
		if (i == Feat2mat.end()){
			Feat2mat[yy] = (int)mat.size();
			mat.push_back(empty);
			FeatureNs.push_back(yy);
			i = Feat2mat.find(yy);
			//cerr << "Could not find entry " << yy << " in registered subset\nAborting.";
			//std::exit(23);
		}
		if ((*i).second > (int)mat.size()){
			cerr << "implied row index larger than high level mat!\nAborting.."; std::exit(25);
		}
		mat[(*i).second][j] += v/div;
	}
}
void HMat::print(ofstream& O){
	O << LvlName << "\t";
	for (unsigned int i = 0; i < SampleNs.size(); i++){
		O << "\t" << SampleNs[i];
	}
	for (unsigned int i = 0; i < FeatureNs.size(); i++){
		O << "\n" << FeatureNs[i];
		for (unsigned int j = 0; j < SampleNs.size(); j++){
			O << "\t" << mat[i][j];
		}
	}
}



VecFiles::VecFiles(const string inF, const string outF, const string xtra) :
infiles(0){
}

int VecFiles::getIdx(const string&){
	int idx(-1);
	return idx;
}


void VecFiles::readVecFile(const string inF){
	string line;
	ifstream in(inF.c_str());

	//int ini_ColPerRow(0),cnt(0);

	//string rowID="";
	while(getline(in,line,'\n')) {
		if(line.substr(0,1) == "#" || line.length()<2){continue;}
		string segments;
		//int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		int cnt2(-1);
		int CurIdx(-1);
		while (getline(ss,segments,'\t')) {
			cnt2++;
			if (cnt2==-1){
				//rowID = segments;
				CurIdx=this->getIdx(segments);
				continue;
			}
			mat_fl tmp =  atof(segments.c_str());
			if (tmp==0){continue;}

		}
	}

}
