struct options
{
public:
	options(int argc, char** argv);
	//~options();

	//vars
  string input = "";
  string output = "";
  string mode  = "";
  uint depth = 0;
  uint repeats = 10;
  uint write = 0;
  uint threads = 1;
  bool writeSwap = true;
  bool verbose = false;

};
