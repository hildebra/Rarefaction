struct rareStruct{
	DivEsts* div;
	string cntsName;
	//std::vector<vector<uint>> cnts;
	vector< map< uint, uint>> cnts;
	string skippedNames;
	vector<string> IDs;
};

void binaryStoreSample(vector< vector< string > > & tmpMatFiles, rareStruct* tmpRS, vector<string>& rowNames, string outF, vector<string>& cntsNames, bool reshapeMap = false);
void memoryStoreSample(rareStruct* tmpRS, vector< vector< map<uint, uint > > >& MaRare, vector<string>& cntsNames,  bool reshapeMap = false);

void printRarefactionMatrix(vector< vector< string > > & tmpMatFiles, string outF, int rareDep, vector<string>& cntsNames, vector<string>& rowNames);
void printRarefactionMatrix(vector < vector< map < uint, uint>>>& MaRare, string outF, int rareDep, vector<string>& cntsNames, vector<string>& rowNames);
