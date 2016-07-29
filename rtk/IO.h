#pragma once
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
//#include <iterator>
#include <string>
#include <cstring>
#include <map>
#include <list>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <random>
#include <assert.h>
#include <unordered_map>
#include <future>
#include <mutex>

//#include <tchar.h>
//#include <string.h>

//multi threading
//#include <thread>
//#include <future>


const bool verbose=0;

#define LINELENGTH 20000;

typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
               // e.g. keep one global instance (per thread)
typedef double mat_fl;
typedef float smat_fl;


using namespace std;
typedef unsigned int uint;
typedef unsigned long ulong;

ulong thr_rng(unsigned long,MyRNG&);
std::istream& safeGetline(std::istream& is, std::string& t);
template<typename T> T getMedian(vector<T>& in){
	sort(in.begin(), in.end());
	size_t size = in.size();
	if (size == 0){ return (T)0; }
	if (size == 1){ return (in[0]) ; }
	if (size == 2){ return (in[0] + in[1]) / 2; }
	T median(in[size / 2]);
	if (size % 2 == 0)	{
		median = (in[size / 2 - 1] + in[size / 2]) / 2;
	}
	return median;
}
void lineCntOut(const string inF, const string outF, const string arg4);

inline std::string stringify(double x)
 {
   std::ostringstream o;
   o << x;
   return o.str();
 }
inline std::string itos(int number) {
	std::stringstream ss;
	ss << number;
	return ss.str();
}

class DivEsts{
public:
	DivEsts():richness(0),shannon(0),
		simpson(0),invsimpson(0),chao1(0),eve(0){}
	~DivEsts(){}
	void print2file(const string);
	//data vectors
	vector<long> richness;
	vector<double> shannon,simpson,invsimpson,chao1,eve;
	string SampleName;
};
void printDivMat(const string outF, vector<DivEsts*>&, bool);
void printRareMat(const string outF, vector< map< uint, uint >>& rMat, vector< string >& sampleNames, vector < string >& rowId);
string printSimpleMap(map<uint, uint> vec, string outF, string id, vector<string> rowNames);
void reassembleTmpMat(vector<string> inF, vector< string > rowNames,vector< string > colNames, string outF);

class smplVec{
public:
	smplVec(const string, const int);
	smplVec(const vector<mat_fl>&, const int);
	~smplVec(){
		//delete[] arr;
	}
	void rarefy(long,string o,int rep,DivEsts*, vector<map<uint, uint>>& RareSample,
		string& retCntsSampleName, string& skippedSample, vector<vector<uint>>* , int=0,bool=false, bool=false);
	long getRichness(const vector<unsigned int>& cnts);
	//int maxSiz(){return vector<unsigned short>::max_size();}
	vector < string > getRowNames(){ return(IDs); }

private:
	int binarySearch(vector<float>,const float x);
	//void shuffle();
	void shuffle_singl();

	//diversity indices
	//method: 1=shannon, 2=simpson, 3=invsimpson
	vector <double> calc_div(const vector<uint>& , int meth, float base=2.718282f);
	double calc_chao1(const vector<uint>& , int corrBias); //corrBias: 0/1
	double calc_eveness(const vector<uint>& );

	void print2File(const vector<unsigned int>&,const string);
	//unsigned short * arr;
	vector<string> IDs;
	vector<unsigned int> arr;
	double totSum;
	vector<MyRNG> rng_P;
	MyRNG rng;
	int num_threads;
	long richness;
	double Shannon;
	int numFeatures;

	//vector<float> vec;
};

vector<mat_fl> computeChao2(vector<vector<uint>>& abundInRow);
void computeICE(vector<vector<uint>>& abundInRow);
void writeChao2(vector<mat_fl>&, string );
