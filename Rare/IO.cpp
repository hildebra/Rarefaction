#include "IO.h"


void lineCntOut(const string inF, const string outF, const string arg4){
	ifstream in(inF.c_str());
	ofstream out(outF.c_str(), ios::out);
	if (!in){ cerr << "Can't open infile " << inF << endl; std::exit(99); }
	if (!out){ cerr << "Can't open outfile " << outF << endl; std::exit(99); }
	//read file that contains nicely ordered all indexes of lines to be extracted
	string line;
	vector<uint> srtTar;
	ifstream idxS(arg4.c_str());
	if (!idxS){ cerr << "Can't open outfile " << arg4 << endl; std::exit(99); }
	while (getline(idxS, line, '\n')) {
		if (line[0] == '>'){
			line.erase(0,1);
		}
		srtTar.push_back(stoi(line));
	}
	idxS.close();
	//sort ascending
	sort(srtTar.begin(), srtTar.end());

	//sort through new file
	if (!out){ cerr << "Can't open outfile " << outF << endl; std::exit(99); }
	int cnt(1); uint j(0);
	while (getline(in, line, '\n')) {
		if (cnt == srtTar[j]){
			out << line << endl;
			uint cur = srtTar[j];
			while (srtTar[j] == cur){ j++; }
			if (j == srtTar.size()){ break; }
		}
		cnt++;
	}

	in.close(); out.close();
	if (j != srtTar.size()){
		cerr << "Missed " << (srtTar.size() - j) << " entries." << endl;
	}
}
//****************************  smplVec::smplVec ***********
smplVec::smplVec(const vector<mat_fl>& vec, const int nt) :IDs(0),totSum(0),
num_threads(nt), richness(-1), Shannon(-1.f){
	double cumSum(0.f);
	for (uint i = 0; i < vec.size(); i++){
		cumSum += vec[i];
	}
	if (verbose){ cerr << (long)cumSum << " allocating "; }
	//arr = (int*) malloc((int)cumSum * sizeof(int));
	//arr = new unsigned short[(int)cumSum];
	arr.resize((long)cumSum);
	if (verbose){ cerr << "memory"; }
	totSum = cumSum;
	long k(0); uint posInVec(-1);
	//numFeatures = 0;
	for (size_t i = 0; i< vec.size(); i++){
		//if (vec.size()-i<10000){cerr<<i<<" ";}
		long maxG = (long)vec[i];

		posInVec++;
		if (maxG == 0){ continue; }//not really a feature, doesnt need ot be counted as cat

		maxG += k; //bring up to current level
		for (; k<maxG; k++){
			arr[k] = posInVec;
		}
		//numFeatures++;

	}
	posInVec++;
	numFeatures = posInVec;
	if (verbose){ cerr << "..\n"; }
}
smplVec::smplVec(const string inF, const int nt) :IDs(0),totSum(0), num_threads(nt),
	richness(-1),Shannon(-1.f) {

	vector<double> vec;
	ifstream in(inF.c_str());
	string line; double cumSum(0.f);
	while(getline(in,line,'\n')) {
		string ID; float num;
		stringstream ss(line);
		ss>>ID; ss>>num;
		cumSum += num;
		vec.push_back(num); IDs.push_back(ID);
	}
	in.close();
	//cerr<<"tt";std::vector<unsigned short> v((int)cumSum);
	if (verbose){cerr<<(long)cumSum<<" allocating ";}
	//arr = (int*) malloc((int)cumSum * sizeof(int));
	//arr = new unsigned short[(int)cumSum];
	arr.resize((long)cumSum);
	if (verbose){cerr<<"memory";}
	totSum = cumSum;
	long k(0); uint posInVec(0);
	for (size_t i = 0; i< vec.size(); i++){
		//if (vec.size()-i<10000){cerr<<i<<" ";}
		long maxG = (long)vec[i];
		maxG += k;
		if (maxG == 0){ continue; }//not really a feature, doesnt need ot be counted as cat
		for (; k<maxG; k++){
			arr[k] = posInVec;
		}
		posInVec++;
	}
	numFeatures = posInVec;
	if (verbose){cerr<<"..\n";}
}

void smplVec::rarefy(long dep, string ofile, int rep,
					DivEsts* divs, std::vector<vector<uint>> & RareSample,
					string& retCntsSampleName, string& skippedSample,
					int writes,bool write, bool fillret){
	if (dep>totSum){
		skippedSample = divs->SampleName;
		if (verbose){cout<<"skipped sample, because rowSums < depth \n";}
		return;
	}
	long curIdx=(long)totSum+1;

	for (int curRep=0;curRep<rep;curRep++){
		if(curIdx+dep >= (long) totSum){
			if (verbose){cerr<<"shuffle \n";}		shuffle_singl();		if (verbose){cerr<<"shed\n";}
			curIdx=0;
		}

		//count up
		vector<unsigned int> cnts(numFeatures, 0);
		for (long i=(0+curIdx);i<(dep+curIdx);i++){
			cnts[arr[i]]++;
		}
		curIdx += dep;
		string t_out = ofile;
		if (rep!=1){
			std::ostringstream oss;
			oss<<curRep;
			t_out += "_" +oss.str();
		}
		if (curRep < writes && write){
			print2File(cnts,t_out);
		}
		if (curRep < writes && fillret) {
			RareSample.push_back(cnts);

			if(curRep == 0){
				retCntsSampleName = divs->SampleName; // safe the sample name as well
			}
		}
		richness = 0;
		divs->richness.push_back(this->getRichness(cnts));
		vector<double> three = this->calc_div(cnts, 4);
		divs->shannon.push_back(three[0]);
		divs->simpson.push_back(three[1]);
		divs->invsimpson.push_back(three[2]);
		divs->chao1.push_back(this->calc_chao1(cnts,1));
		divs->eve.push_back(this->calc_eveness(cnts));
		richness = 0;
	}
}

long smplVec::getRichness(const vector<unsigned int>& cnts){
	for (size_t i = 0; i<cnts.size(); i++){
		//out<<IDs[i]<<"\t"<<cnts[i]<<endl;
		if (cnts[i]>0){
			richness++;
		}
	}
	return richness;
}


void smplVec::print2File(const vector<unsigned int>& cnts,const string t_out){
	richness=0;
	ofstream out(t_out.c_str());
	for (size_t i=0;i<cnts.size();i++){
		//out<<IDs[i]<<"\t"<<cnts[i]<<endl;
		if (cnts[i]>0){
			richness++;
			out<<IDs[i]<<"\t"<<cnts[i]<<endl;
		}
	}
	out.close();
	cout<<"Richness: "<<richness<<endl;
	//return richness;
}

ulong thr_rng(unsigned long pos,MyRNG& rng) {
    std::uniform_int_distribution<unsigned long> uint_distx(0,pos);
	return uint_distx(rng);
}

/*
void smplVec::shuffle(){
	time_t seed_val=time(NULL);           // populate somehow
	rng_P.resize(num_threads);
	for (long t=0; t<num_threads; t++){ //initialize N random number generators
		rng_P[t].seed((long)seed_val-(t*13));
	}
	std::vector<std::future<ulong>> futures(num_threads);
	//auto *thr = new async[num_threads];
	//vector<unsigned long> rngs(num_threads,0);
	ulong trng;

	unsigned long pos = (unsigned long)totSum- 1;
	while ( pos > 0) {
		//cout << pos<<" ";
		        //Launch a group of threads
        for (int t = 0; t < num_threads; ++t) {
            futures[t] = async(thr_rng,pos-t,rng_P[t]);
        }
		for (int t = 0; t < num_threads; ++t) {
			trng = futures[t].get();
			unsigned int temp = arr[pos] ;
			arr[pos] = arr[trng];
			arr[trng] = temp;
			pos--;
        }
		//std::uniform_int_distribution<unsigned long> uint_distx(0,pos);
		//unsigned long j = uint_distx(rng);
		//swap(arr[i],arr[j]);
	}
	//delete [] thr;
	if (verbose){cerr<<"fini";}
}
*/
void smplVec::shuffle_singl(){
	time_t seed_val=time(NULL);           // populate somehow
	rng.seed((long)seed_val);

	for (unsigned long i = (unsigned long)totSum- 1; i > 0; i--) {
		std::uniform_int_distribution<unsigned long> uint_distx(0,i);
		unsigned long j = uint_distx(rng);
		unsigned int temp = arr[i] ;
		arr[i] = arr[j];
		arr[j] = temp;
		//swap(arr[i],arr[j]);
	}
	if (verbose){cerr<<"fini";}
}



int smplVec::binarySearch( vector<float> vec, const float toFind)
{
	int len =(int) vec.size();
    // Returns index of toFind in sortedArray, or -1 if not found
    int low = 0;
    int high = len - 1;
    int mid;

    float l = vec[low];
    float h = vec[high];

    while (l < toFind && h > toFind) {
        mid = (low + high)>>1;

        float m = vec[mid];

        if (m < toFind) {
			if (vec[mid+1] > toFind){return mid;}
            l = vec[low = mid + 1];
        } else if (m > toFind) {
			if (vec[mid-1] < toFind){return mid-1;}
            h = vec[high = mid - 1];
        } else {
            return mid;
        }
		if (mid==len || mid==0){return mid;}
    }
    if (vec[low] == toFind)
        return low;
    else
        return -1; // Not found
}

/*
	// my own implementation of Chao2
	// http://viceroy.eeb.uconn.edu/EstimateSPages/EstSUsersGuide/EstimateSUsersGuide.htm#AppendixB
	#for incidence data
	.chao2 = function(M,bias.corr=T,conf.int=F){
		#browser()
		m = dim(M)[2]
		M[M>1]=1; #convert to incidence data
		pool = apply(M,1,sum)
		s_obs = apply(M>0,2,sum)
		singlIdx = pool==1
		doublIdx = pool==2
		Q1 = apply(M[singlIdx,],2,sum)
		Q2 = apply(M[doublIdx,],2,sum)
		if (!bias.corr){
			est = s_obs + (Q1^2)/(2*Q2)
		} else {
			est = s_obs + ((m-1)/m) *(Q1*(Q1-1))/(2*(Q2+1));
		}
		if (conf.int){
			stop("TODO")
			M="?"
			P = exp(-M/Sobs)
			P1 =Sobs/(1-P)
			P2 = 1.96*sqrt((Sobs*P)/(1-P))
			low =  P1 - P2
			idx = low<Sobs
			low[idx] = Sobs[idx]
			hi = P1 + P2
		}

		return (est)
	}
	/*/

double smplVec::calc_chao1(const vector<uint> & vec,int corrBias=1){
	double Sobs((double)richness);
	double singl(0); double doubl(0);
	for (size_t i=0;i<vec.size();i++){
		if (vec[i]==1){singl++;}
		else if (vec[i]==2){doubl++;}
	}
	double est=0.0;
	if (corrBias==0){
		est = float( Sobs + (singl*singl)/(2*doubl) );
	} else {
		est = float( Sobs + (singl*(singl-1))/(2*(doubl+1)) );
	}
	/*if (conf.int){
	N = apply(M,2,sum)
	P = exp(-N/Sobs)
	P1 =Sobs/(1-P)
	P2 = 1.96*sqrt((Sobs*P)/(1-P))
	low =  P1 - P2
	idx = low<Sobs
	low[idx] = Sobs[idx]
	hi = P1 + P2
	}*/
	return (est);
}


vector<double> smplVec::calc_div(const vector<uint>& vec,int meth=1, float base){
	double sum = 0;
	for (size_t i=0; i<vec.size();i++){sum+=(double)vec[i];}
	vector<double> x(vec.begin(),vec.end());
	for (size_t i=0; i<x.size();i++){x[i] /= sum;}
	bool doexp = false;
	if (base <= 2.718284f && base >= 2.718280f){ // account for machine imprecission
		doexp = true;
	}
	vector<double> H(3, 0.0); double H1(0.0), H2(0.0), H3(0.0);
	if (meth == 1 || meth == 4){
		if (doexp){
			for (size_t i=0; i<x.size();i++){if (x[i]>0){H1 += x[i] * -log(x[i])  ;}}
		} else {
			float div = -log10(base);
			for (size_t i = 0; i<x.size(); i++){ if (x[i]>0){ H1 += x[i] * log10(x[i]) / div; } }
		}
		Shannon = H1;
	}
	if (meth == 3 || meth == 4 || meth == 2) {
		for (size_t i = 0; i<x.size(); i++){ H2 += x[i] * x[i]; }
		H3 = H2;
		H2 = 1 - H2;//simpson
		H3 = 1 / H3; //invsimpson
	}
	//for (size_t i=0; i<x.size();i++){H += x[i];}
	//if (meth == (int)2) {		H = 1 - H;	}else if (meth == 3)		H = 1/H;	}
	H[0] = H1; H[1] = H2; H[2] = H3;

	return(H);
}
double smplVec::calc_eveness(const vector<uint>& vec){
	//double sha = calc_div(vec,1);
	if (Shannon == -1.f){ vector<double> tm = calc_div(vec, 1); }
	return(Shannon / log((double)richness));
}



void DivEsts::print2file(const string file){
	if (richness.size()<1){return;}
	ofstream out(file.c_str());
	if (!out){ cerr << "Couldn't open diversity estimate file " << file << endl; std::exit(99); }
	out<<"Richness\t"<<richness[0];
	for (size_t i=1; i<richness.size();i++){
		out << "\t"<<richness[i];
	}
	out<<"\nShannon\t"<<shannon[0];
	for (size_t i=1; i<shannon.size();i++){
		out << "\t"<<shannon[i];
	}
	out<<"\nSimpson\t"<<simpson[0];
	for (size_t i=1; i<simpson.size();i++){
		out << "\t"<<simpson[i];
	}
	out<<"\nInv. Simpson\t"<<invsimpson[0];
	for (size_t i=1; i<invsimpson.size();i++){
		out << "\t"<<invsimpson[i];
	}
	out<<"\nChao1\t"<<chao1[0];
	for (size_t i=1; i<chao1.size();i++){
		out << "\t"<<chao1[i];
	}
	out<<"\nEveness\t"<<eve[0];
	for (size_t i=1; i<eve.size();i++){
		out << "\t"<<eve[i];
	}
	out.close();
}
void printDivMat(const string outF, vector<DivEsts*>& inD){
	ofstream out(outF.c_str());
	if (!out){ cerr << "Couldn't open diversity estimate matrix " << outF << endl; std::exit(99); }
	out << "Smpl\tRichness\tShannon\tSimpson\tInv. Simpson\tChao1\tEveness\n";
	for (size_t i = 0; i < inD.size(); i++){
		if (inD[i] == NULL){ cerr << "Empty vector at index " << i << "in div mat building.\n";
			out << "-1\t-1\t-1\t-1\t-1\t-1\n";
			continue;
		}
		out << inD[i]->SampleName << "\t";
		out << getMedian(inD[i]->richness) << "\t";
		out << getMedian(inD[i]->shannon) << "\t";
		out << getMedian(inD[i]->simpson) << "\t";
		out << getMedian(inD[i]->invsimpson) << "\t";
		out << getMedian(inD[i]->chao1) << "\t";
		out << getMedian(inD[i]->eve) << "\n";
	}
	out.close();
}


std::istream& safeGetline(std::istream& is, std::string& t)
{
	t.clear();
	//from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
	// The characters in the stream are read one-by-one using a std::streambuf.
	// That is faster than reading them one-by-one using the std::istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.

	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();

	for (;;) {
		int c = sb->sbumpc();
		switch (c) {
		case '\n':
			return is;
		case '\r':
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:
			// Also handle the case when the last line has no line ending
			if (t.empty())
				is.setstate(std::ios::eofbit);
			return is;
		default:
			t += (char)c;
		}
	}
}
