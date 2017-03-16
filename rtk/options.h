#pragma once
//#include "IO.h"

struct options
{
public:
	options(int argc, char** argv);
	options(std::string, std::string , int repeats, std::vector<long> depth, int NoOfMatrices, bool verbose, unsigned int threads);
	void print_rare_details();
	//~options();

	//vars
  std::string input = "";
  std::string output = "";
  std::string mode  = "";
  std::string referenceDir = "";
  std::string referenceFile = "";
  std::string map = "";
  std::vector<long> depth;
  long depthMin;
  uint repeats = 10;
  uint write = 0;
  uint threads = 1;
  bool writeSwap = true;
  bool verbose = false;

  std::string modDB;
  int modRedund;
  float modEnzCompl;
  float modModCompl;
  bool modWrXtraInfo;
  bool modCollapse;
  bool calcCoverage;

  std::string xtra;
};
