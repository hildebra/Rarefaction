#pragma once
#include "ClStr2Mat.h"
#include "options.h"

struct rareStruct{
	int i;
	DivEsts* div;
	string cntsName;
	vector<vector<rare_map>> cnts;
	string skippedNames;
	vector<string> IDs;
	
};

struct job {
  std::future <rareStruct*> fut;
  bool inUse = false;
};

void binaryStoreSample(vector<vector< vector< string >> > & , rareStruct* , 
	vector<string>& , string , vector<string>& , bool reshapeMap = false);
void memoryStoreSample(rareStruct* tmpRS, vector< vector< vector< rare_map >> >& MaRare,  vector<string>& cntsNames, bool reshapeMap);

int rarefyMain(options* opts,  string mode,
	vector<vector<mat_fl>> rmatrix,
	vector< string > cnames , vector< string > rnames ,
	vector<DivEsts*>&  divvs,
	std::vector<vector<vector<rare_map>>> &MaRare,
	std::vector<string>& cntsNames,
	std::vector<string>& skippedNames,
	std::vector<vector<mat_fl>>& ACE,
	std::vector<vector<mat_fl>>& ICE,
	std::vector<vector<mat_fl>>& chao2,
	std::vector<string>& rowNames,
	bool transpose);
