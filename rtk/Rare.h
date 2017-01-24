#pragma once
#include "ClStr2Mat.h"
#include "options.h"

struct rareStruct{
	int i;
	DivEsts* div;
	string cntsName;
	//std::vector<vector<uint>> cnts;
	vector<rare_map> cnts;
	string skippedNames;
	vector<string> IDs;
	
};

struct job {
  std::future <rareStruct*> fut;
  bool inUse = false;
};

void binaryStoreSample(vector< vector< string > > & tmpMatFiles, rareStruct* tmpRS, 
	vector<string>& rowNames, string outF, vector<string>& cntsNames, bool reshapeMap = false);
void memoryStoreSample(rareStruct* tmpRS, vector< vector< rare_map > >& MaRare, 
	vector<string>& cntsNames,  bool reshapeMap = false);

void printRarefactionMatrix(vector< vector< string > > & tmpMatFiles, string outF, 
	int rareDep, vector<string>& cntsNames, vector<string>& rowNames);
void printRarefactionMatrix(const vector < vector< rare_map >>& MaRare, string outF,
	int rareDep, vector<string>& cntsNames, vector<string>& rowNames);
