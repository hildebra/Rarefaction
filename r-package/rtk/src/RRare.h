#pragma once
#include "ClStr2Mat.h"
#include "options.h"

int rarefyMain(string inF, string outF, string mode ,
    int repeats, long rareDep, unsigned int numThr, bool verbose,
    vector<vector<mat_fl>> rmatrix,
    vector< string > cnames , vector< string > rnames, vector<DivEsts*> * divvs,
    std::vector<vector<rare_map>> &retCnts,
  	std::vector<string>& retCntsSampleNames,
  	std::vector<string>& skippedSamples,
    std::vector<double>& ACE, std::vector<double>& ICE, std::vector<double>& chao2,
  	std::vector<string>& rowNames,
  	int NoOfMatrices, bool transpose);
