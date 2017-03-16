#include "RRare.h"

const char* rar_ver="0.93r";
/*
struct cDR{
	DivEsts* div;
	std::vector<vector<vector<uint>>> retCnts;
	string retCntsSampleName;
	vector<rare_map> RareSample;
	string skippedSample;
};*/


rareStruct* calcDivRar(int i, Matrix* Mo, DivEsts* div, options* opts,
        vector<vector<vector<uint>>>* abundInRow, vector<vector<vector<uint>>>* occuencesInRow){

    smplVec* cur        = Mo->getSampleVec(i);
    string curS         = Mo->getSampleName(i);
    div->SampleName     = curS;
    std::vector<vector<uint>> cnts;
    vector<vector<rare_map>> cntsMap(opts->depth.size());
    vector<string> cntsName(opts->depth.size());
    string skippedNames;
    bool wrAtAll(opts->write > 0);

    cur->rarefy(opts->depth, opts->output, opts->repeats,
            div, cntsMap, cntsName, skippedNames, abundInRow, occuencesInRow,
            opts->write, false, wrAtAll);

    // put everything in our nice return container
    rareStruct* tmpRS       = new rareStruct();
    tmpRS->div              = div;
    tmpRS->cnts             = cntsMap;
    tmpRS->cntsName         = cntsName;
    tmpRS->skippedNames     = skippedNames;
    tmpRS->i                = i;

    delete cur;
    return tmpRS;
}


rareStruct* calcDivRarVec(int i, vector<string> fileNames, DivEsts* div, options* opts,
        vector<vector<vector<uint>>>* abundInRow, vector<vector <vector<uint>>>* occuencesInRow){
    smplVec* cur = new smplVec(fileNames[i],4);

    std::vector<vector<uint>> cnts;
    vector<vector<rare_map>> cntsMap(opts->depth.size());
    vector<string> cntsName(opts->depth.size());
    string skippedNames;
    bool wrAtAll(opts->write > 0);
    cur->rarefy(opts->depth, opts->output, opts->repeats,
            div, cntsMap, cntsName, skippedNames, abundInRow, occuencesInRow,
            opts->write, false, wrAtAll);

    rareStruct* tmpRS       = new rareStruct();
    tmpRS->div              = div;
    tmpRS->cnts             = cntsMap;
    tmpRS->cntsName         = cntsName;
    tmpRS->skippedNames     = skippedNames;
    tmpRS->IDs              = cur->getRowNames();
    tmpRS->i                = i;

    delete cur;
    if( remove( fileNames[i].c_str() ) != 0 ){
        //cerr << "Swap Mode: Error deleting file: " << fileNames[i] << std::endl;
    }
    return tmpRS;
}


void rareLowMem(options* opts,     
            vector<DivEsts*>&  divvs,
            std::vector<vector<vector<rare_map>>> &MaRare,
            std::vector<vector<string>>& cntsNames,
            std::vector<string>& skippedNames,
            std::vector<vector<mat_fl>>& ACE,
            std::vector<vector<mat_fl>>& ICE,
            std::vector<vector<mat_fl>>& chao2,
            std::vector<string>& rowNames,
            bool transpose){
                
    vector<string> fileNames;
    Matrix* Mo 	= new Matrix(opts->input, opts->output, "", fileNames, false, true);
    vector < string > SampleNames 	= Mo->getSampleNames();
    rowNames 		= Mo->getRowNames();
    
    // abundance vectors to hold the number of occurences of genes per row
    // this will be used for ICE, ACE and or Chao2 estimation
    vector<vector<vector<uint>>> occuencesInRow(opts->depth.size(), vector<vector<uint>>(opts->repeats, vector<uint>(Mo->rowNum(),0)));
    vector<vector<vector<uint>>> abundInRow(opts->depth.size(), vector<vector<uint>>(opts->repeats, vector<uint>(Mo->rowNum(),0)));

    for(uint i = 0; i < opts->depth.size(); i++){
            if (opts->depth[i] < 1.) {
                // rarefy to smallest colSum
                opts->depth[i] = (uint)round(opts->depth[i] * Mo->getMinColSum());
                if (opts->depth[i] == 0.0) {
                    //cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
                    return;
                }
            } 
        }
    size_t smpls = Mo->smplNum();
    divvs.resize(smpls, NULL);
    delete Mo;
    
   

    vector < vector < vector < string >> > tmpMatFiles(opts->depth.size(), vector<vector <string>>(opts->write));

    vector < job > slots(opts->threads);

// now start a async in each slot
    uint i          = 0; 
    bool breakpoint(true);

    while (breakpoint) {
        // check for any finished jobs
        for( uint j = 0; j < slots.size(); j++ ){
            if( i >= smpls){
                breakpoint = false;
                // break in case we have more slots than work
                break;
            }

            if(slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(20)) == std::future_status::ready){
                // move the information
                rareStruct* tmpRS;
                tmpRS               = slots[j].fut.get();
                divvs[tmpRS->i]     = tmpRS->div;
                string curS 	    = SampleNames[tmpRS->i];
                // add the matrices to the container
                if (opts->write > 0) {
                    if (opts->writeSwap) {
                        memoryStoreSample(opts, tmpRS, MaRare, cntsNames, true);
                    }
                    else {
                        memoryStoreSample(opts, tmpRS, MaRare, cntsNames, true);
                    }
                }

                delete tmpRS;
                // free slot
                slots[j].inUse = false;
            }

            // open new slots
            if( slots[j].inUse == false){
                slots[j].inUse = true;
                // launch an async task
                DivEsts * div   = new DivEsts();
                div->SampleName = SampleNames[i];
                slots[j].fut    = async(std::launch::async, calcDivRarVec, i, fileNames, div, opts, &abundInRow, &occuencesInRow);
                

                i++;

            }

        }


    }

    // wait for the threads to finish up.
    for(uint j = 0; j < slots.size(); j++){
        if(slots[j].inUse == false ){
            // only copy if there is work to be done
            continue;
        }
        slots[j].fut.wait();
        // move the information
        rareStruct* tmpRS;
        tmpRS               = slots[j].fut.get();
        divvs[tmpRS->i]     = tmpRS->div;
        string curS 	    = SampleNames[tmpRS->i];

        // add the matrices to the container
        if (opts->write > 0) {
            if (opts->writeSwap) {
                memoryStoreSample(opts, tmpRS, MaRare, cntsNames, true);
            }
            else {
                memoryStoreSample(opts, tmpRS, MaRare, cntsNames, true);
            }
        }

        delete tmpRS;
        // free slot
        slots[j].inUse = false;
    }



    // compute chao2, ACE, ICE and write to file
    computeChao2(chao2, abundInRow);
    computeCE(ICE, abundInRow);
    computeCE(ACE, occuencesInRow);


}


void binaryStoreSample(options* opts, vector<vector<vector< string >>>& tmpMatFiles, rareStruct* tmpRS, vector<string>& rowNames, string outF,  vector<vector<string>>& cntsNames, bool reshapeMap){
    // store vectors of rarefied matrix on hard disk for memory reduction
    if(reshapeMap){
        vector < string > rowIDs = tmpRS->IDs;
        vector < uint > nrowIDs(rowIDs.size());
        // convert ids into integer vector, for remapping the values
        for(uint i = 0; i < rowIDs.size(); i++){
            nrowIDs[i] = std::stoi(rowIDs[i]);
        }
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            // reshape each vector, as some are zero, and we need to rematch values and rows
            // we use the row Ids which we created correctly when splitting the vector from the input file
            rare_map tmpVec;
            for (auto const& x : tmpRS->cnts[i][ii]){
                tmpVec[nrowIDs[x.first]] = x.second;
            }
            string vecLocation = printSimpleMap(tmpVec,	outF + "tmp_"  + std::to_string(ii)+"_" + std::to_string(i) + tmpRS->cntsName[i] + ".binary",	tmpRS->cntsName[i], rowNames);

            tmpMatFiles[i][ii].push_back(vecLocation);
        }
        }
    }else{
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            string vecLocation = printSimpleMap(tmpRS->cnts[i][ii],	outF + "tmp_" + std::to_string(ii) + "_" + std::to_string(i) + tmpRS->cntsName[i] + ".binary",	tmpRS->cntsName[i], rowNames);
            tmpMatFiles[i][ii].push_back(vecLocation);
        } }
    }
    // save sample name
    for(uint di = 0; di < opts->depth.size(); di++){
        if(tmpRS->cntsName[di].size() != 0){
            cntsNames[di].push_back(tmpRS->cntsName[di]);

        }
    }
}

void memoryStoreSample(options* opts, rareStruct* tmpRS, vector< vector< vector< rare_map >> >& MaRare,  vector<vector<string>>& cntsNames, bool reshapeMap){
    if(reshapeMap){
        vector < string > rowIDs = tmpRS->IDs;
        vector < uint > nrowIDs(rowIDs.size());
        // convert ids into integer vector, for remapping the values
        for(uint i = 0; i < rowIDs.size(); i++){
            nrowIDs[i] = std::stoi(rowIDs[i]);
        }
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            // reshape each vector, as some are zero, and we need to rematch values and rows
            // we use the row Ids which we created correctly when splitting the vector from the input file
            rare_map tmpVec;
            for (auto const& x : tmpRS->cnts[i][ii]){
                tmpVec[nrowIDs[x.first]] = x.second;
            }
            MaRare[i][ii].push_back(tmpVec);
        }
    }
    }else{
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            MaRare[i][ii].push_back(tmpRS->cnts[i][ii]);
            }
        }
    }
    // save sample name
    for(uint di = 0; di < opts->depth.size(); di++){
        if(tmpRS->cntsName[di].size() != 0){
            cntsNames[di].push_back(tmpRS->cntsName[di]);
        }
    }
}










//int main(int argc, char* argv[])
int rarefyMain(options* opts,  string mode,
	vector<vector<mat_fl>> rmatrix,
	vector< string > cnames , vector< string > rnames ,
	vector<DivEsts*>&  divvs,
	std::vector<vector<vector<rare_map>>> &MaRare,
	std::vector<vector<string>>& cntsNames,
	std::vector<string>& skippedNames,
	std::vector<vector<mat_fl>>& ACE,
	std::vector<vector<mat_fl>>& ICE,
	std::vector<vector<mat_fl>>& chao2,
	std::vector<string>& rowNames,
	bool transpose)
{
	// compatibility to main rare software
	bool writeFiles = false;

	MyRNG rng;


	if (mode == "memory"){
		Matrix* Mo;
		if(opts->input != ""){
			  Mo = new Matrix(opts->input, "");//no arg for outfile &  hierachy | gene subset
		}else{
		  Mo = new Matrix(); // empty matrix for filling later
		  int nr = rmatrix.size();
		  for(int i=0; i < nr; i++){
		    Mo->addRow(rmatrix[i]);
		  }
		  Mo->setSampleNames(cnames);
		  Mo->setRowNames(rnames);
		}

		if(transpose == true){
			Mo->transpose();
		}

        // transform all percentage sizes into correct values
        for(uint i = 0; i < opts->depth.size(); i++){
            if (opts->depth[i] < 1.) {
                // rarefy to smallest colSum
                opts->depth[i] = (uint)round(opts->depth[i] * Mo->getMinColSum());
                if (opts->depth[i] == 0.0) {
                    return 0; // without warning return, sad but should never happen
                    // as we check this in R as a sanity check
                }
            } 
        }
		rowNames = Mo->getRowNames();
		// abundance vectors to hold the number of occurences of genes per row
		// this will be used for Chao2 estimation
        vector<vector<vector<uint>>> occuencesInRow(opts->depth.size(), vector<vector<uint>>(opts->repeats, vector<uint>(Mo->rowNum(),0)));
        vector<vector<vector<uint>>> abundInRow(opts->depth.size(), vector<vector<uint>>(opts->repeats, vector<uint>(Mo->rowNum(),0)));
            vector < vector < vector < string >> > tmpMatFiles(opts->depth.size(), vector<vector <string>>(opts->write));
        
        divvs.resize(Mo->smplNum(), NULL);
        
        
        // vector keeping all the slots
        vector < job > slots(opts->threads);

        // now start a async in each slot
        uint i          = 0; 
        size_t smpls = Mo->smplNum();
        bool breakpoint(true);

    while (breakpoint) {
        // check for any finished jobs
        for( uint j = 0; j < slots.size(); j++ ){
            if( i >= smpls){
                breakpoint = false;
                // break in case we have more slots than work
                break;
            }

            if(slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(20)) == std::future_status::ready){

                // move the information
                rareStruct* tmpRS;
                tmpRS               = slots[j].fut.get();
                divvs[tmpRS->i]     = tmpRS->div;
                string curS         = Mo->getSampleName(tmpRS->i);

                // add the matrices to the container
                if (opts->write > 0) {
                    if (opts->writeSwap) {
                       memoryStoreSample(opts, tmpRS, MaRare, cntsNames, false);
                    }
                    else {
                        memoryStoreSample(opts, tmpRS, MaRare, cntsNames, false);
                    }
                }

                delete tmpRS;
                // free slot
                slots[j].inUse = false;
            }

            // open new slots
            if( slots[j].inUse == false){

                slots[j].inUse = true;
                // launch an async task
                DivEsts * div   = new DivEsts();
                slots[j].fut    = async(std::launch::async, calcDivRar, i, Mo, div, opts, &abundInRow, &occuencesInRow);

                i++;
            }
        }
    }

    // wait for the threads to finish up.
    for(uint j = 0; j < slots.size(); j++){
        if(slots[j].inUse == false ){
            // only copy if there is work to be done
            continue;
        }
        slots[j].fut.wait();
        // move the information
        rareStruct* tmpRS;
        tmpRS               = slots[j].fut.get();
        divvs[tmpRS->i]     = tmpRS->div;
        string curS         = Mo->getSampleName(tmpRS->i);

        // add the matrices to the container
        if (opts->write > 0) {
            if (opts->writeSwap) {
               memoryStoreSample(opts, tmpRS, MaRare, cntsNames, false);
            }
            else {
                memoryStoreSample(opts, tmpRS, MaRare, cntsNames, false);
            }
        }

        delete tmpRS;
        // free slot
        slots[j].inUse = false;
    }

		delete Mo;

		// compute chao2, ACE, ICE and write to file
		computeChao2(chao2, abundInRow);
		computeCE(ICE, abundInRow);
		computeCE(ACE, occuencesInRow);


	}else if(mode == "swap"){
		rareLowMem(opts, divvs, MaRare, cntsNames,
            skippedNames, ACE, ICE, chao2, rowNames,
            transpose);
	}


	return 0;
}




