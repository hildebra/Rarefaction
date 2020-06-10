#pragma once
//#include "IO.h"

struct options
{
public:
	options(int argc, char** argv);
	options(std::string, std::string , int repeats, std::vector<double> depth, 
            unsigned int seed,
            int NoOfMatrices, 
		bool verbose, unsigned int threads);
	void print_rare_details();
	//~options();

	//vars
  std::string input = "";
  std::string output = "";
  std::string mode  = "";
  std::string referenceDir = "";
  std::string referenceFile = "";
  std::string map = "";
  std::vector<double> depth;
  
  long depthMin;
  unsigned int repeats = 10;
  unsigned int seed = 1;
  unsigned int write = 0;
  unsigned int threads = 1;
  bool writeSwap = true;
  bool verbose = false;
  bool oldMapStyle = true;

    std::string modDB;
    int modRedund;
    float modEnzCompl;
    float modModCompl;
    bool modWrXtraInfo;
    bool modCollapse;
	bool calcCoverage;
	bool calcCovMedian;
	bool check4idxMatch;//assummes tab separated row name, that is number (idx); used in "lineExtr"

  std::string modDescr;
  std::string modHiera;
  std::string xtra;
};
