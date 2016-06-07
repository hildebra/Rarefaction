
#include "ClStr2Mat.h"

int rarefyMain(string inF, string outF, string mode ,
    int repeats, long rareDep, unsigned int numThr, bool verbose,
    bool returnObject, vector<vector<mat_fl>> rmatrix,
    vector< string > cnames , vector< string > rnames , DivEsts * dd, vector<DivEsts*> * divvs,
    std::vector<vector<vector<uint>>> &retCnts,
	std::vector<string>& retCntsSampleNames, std::vector<string>& rowNames,
	int NoOfMatrices, bool transpose);
