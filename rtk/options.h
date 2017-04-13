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
    std::string input;
    std::string output;
    std::string mode;
    std::string referenceDir;
    std::string referenceFile;
    std::string map;
    std::vector<long> depth;
    long depthMin;
    unsigned int repeats ;
    unsigned int write ;
    unsigned int threads ;
    bool writeSwap ;
    bool verbose;

    std::string modDB;
    int modRedund;
    float modEnzCompl;
    float modModCompl;
    bool modWrXtraInfo;
    bool modCollapse;
    bool calcCoverage;

    std::string xtra;
};
