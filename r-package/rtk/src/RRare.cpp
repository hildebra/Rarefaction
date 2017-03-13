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
    string cntsName;
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
    string cntsName;
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

void helpMsg(){
}

/*

void rareLowMem(string inF, string outF, int NoOfMatrices, long arg4, int repeats,
	vector<DivEsts*>&  divvs,
	std::vector<vector<rare_map>>& MaRare,
	std::vector<string>& cntsNames,
	std::vector<string>& skippedNames,
	std::vector<string>& rowNames, int numThr, bool verbose ){
	// this mode takes the file, reads it in memory
	// prints the columns to their own files
	// then it loads those files again and
	// rarefies each column
	// the measures are then combines again.

	//split mat code

	vector<string> fileNames;
	Matrix* Mo 	= new Matrix(inF, outF, "", fileNames, false, true);
	vector < string > SampleNames 	= Mo->getSampleNames();
	rowNames 		= Mo->getRowNames();

	int rareDep 	= arg4;
	if(rareDep == 0){
		// rarefy to smallest colSum
		rareDep = round(0.95 * Mo->getMinColSum());
		//if(rareDep == 0){
		//	cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
		//	return;
		//}
	}
	delete Mo;

	int done = 0; // number of samples processed for multithreading
	uint i = 0;
	std::future<rareStruct*> *tt = new std::future<rareStruct*>[numThr - 1];
	// abundance vectors to hold the number of occurences of genes per row
	// this will be used for Chao2 estimation
	vector<vector<uint>> abundInRow(repeats, vector<uint>(Mo->rowNum(),0));
	vector<vector<uint>> occuencesInRow(repeats, vector<uint>(Mo->rowNum(),0));

	//rarefection code
	while(i < fileNames.size()){

		// allow multithreading
		//if(verbose == true){
		//	cerr << "At Sample " << i+1 << " of " << fileNames.size() << " Samples";
		//}
		uint toWhere = done + numThr - 1;
		if ((uint)((uint)fileNames.size() - 2 ) < toWhere){
			toWhere = fileNames.size() - 2;
		}
		// launch samples in threads
		for (; i < toWhere; i++){
			DivEsts * div 	= new DivEsts();
			div->SampleName = SampleNames[i];
			tt[i - done] = async(std::launch::async, calcDivRarVec, i, fileNames,  div, rareDep,
				 									&abundInRow, &occuencesInRow, outF, repeats, NoOfMatrices);
		}
		// launch one in the mainthread
		DivEsts * div 	= new DivEsts();
		div->SampleName = SampleNames[i];
		rareStruct* tmpRS;
		tmpRS = calcDivRarVec(i, fileNames,  div, rareDep,
													 &abundInRow, &occuencesInRow, outF, repeats, NoOfMatrices);
		i++;

		// process created data, first threads, then main thread
		i = done;
		for (; i < toWhere; i++){
			rareStruct* RSasync;
			RSasync 		= tt[i-done].get();
			divvs.push_back(RSasync->div);
			string curS 	= SampleNames[i];


			// add the matrices to the container
			if(NoOfMatrices > 0){
				for(uint i = 0; i < RSasync->cnts.size(); i++){
					MaRare[i].push_back(RSasync->cnts[i]);
				}
				// save sample name for naming purposes
				if(RSasync->cntsName.size() != 0){
					cntsNames.push_back(RSasync->cntsName);
				}
			}
			delete RSasync;
		}

		// main thread divv push back
		divvs.push_back(tmpRS->div);
		string curS 	= SampleNames[i];

		if(NoOfMatrices > 0){
			for(uint i = 0; i < tmpRS->cnts.size(); i++){
				MaRare[i].push_back(tmpRS->cnts[i]);
			}

			// save sample name for naming purposes
			if(tmpRS->cntsName.size() != 0){
				cntsNames.push_back(tmpRS->cntsName);
			}
		}
		delete tmpRS;
		i++;
		done = i;
	}

}


*/


void binaryStoreSample(vector<vector<vector< string >>>& tmpMatFiles, rareStruct* tmpRS, vector<string>& rowNames, string outF,  vector<string>& cntsNames, bool reshapeMap){
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
            string vecLocation = printSimpleMap(tmpVec,	outF + "tmp_"  + std::to_string(ii)+"_" + std::to_string(i) + tmpRS->cntsName + ".binary",	tmpRS->cntsName, rowNames);
            tmpMatFiles[i][ii].push_back(vecLocation);
        }
        }
    }else{
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            string vecLocation = printSimpleMap(tmpRS->cnts[i][ii],	outF + "tmp_" + std::to_string(ii) + "_" + std::to_string(i) + tmpRS->cntsName + ".binary",	tmpRS->cntsName, rowNames);
            tmpMatFiles[i][ii].push_back(vecLocation);
        } }
    }
    // save sample name
    if(tmpRS->cntsName.size() != 0){
        cntsNames.push_back(tmpRS->cntsName);
    }
}

void memoryStoreSample(rareStruct* tmpRS, vector< vector< vector< rare_map >> >& MaRare,  vector<string>& cntsNames, bool reshapeMap){
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
    if(tmpRS->cntsName.size() != 0){
        cntsNames.push_back(tmpRS->cntsName);
    }
}











//int main(int argc, char* argv[])
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
        
        divvs.resize(Mo->smplNum());

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
                       binaryStoreSample(tmpMatFiles, tmpRS, rowNames, opts->output, cntsNames, false);
                    }
                    else {
                        memoryStoreSample(tmpRS, MaRare, cntsNames, false);
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
               binaryStoreSample(tmpMatFiles, tmpRS, rowNames, opts->output, cntsNames, false);
            }
            else {
                memoryStoreSample(tmpRS, MaRare, cntsNames, false);
            }
        }

        delete tmpRS;
        // free slot
        slots[j].inUse = false;
    }

    // output matrix
    //printDivMat(opts->output, divvs, true, opts);
    for (size_t i = 0; i < divvs.size(); i++) {
        delete divvs[i];
    }

    // write rarefaction matrices to disk
    /*if (opts->write > 0) {
        vector< string > rowNames = Mo->getRowNames();
        if (opts->writeSwap) {
            printRarefactionMatrix(opts, tmpMatFiles, opts->output, cntsNames, rowNames);
        }
        else {
            printRarefactionMatrix(opts, MaRare, opts->output,  cntsNames, rowNames);
        }
    }*/
    



/*
		//Rcout << "threads\n";
		uint i = 0; uint done = 0;
		while ( i < Mo->smplNum()){

			uint toWhere = done+numThr - 1;
			if ((uint)((uint)Mo->smplNum() - 2 ) < toWhere){
					toWhere = Mo->smplNum() - 2;
			}
			for (; i < toWhere; i++){ // with just one thread this is not used?
				DivEsts * div = new DivEsts();
				tt[i - done] = async(std::launch::async, calcDivRar, i, Mo, div, rareDep,
								 	&abundInRow, &occuencesInRow, outF,
									repeats, writeFiles, NoOfMatrices);

			}

			//use main thread to calc one sample as well
			DivEsts * div = new DivEsts();
			
			rareStruct* tmpCDr;
			//tmpCDr = new rareStruct;
			tmpCDr 		= calcDivRar(i, Mo, div, rareDep,  &abundInRow, &occuencesInRow,
									outF, repeats, writeFiles,
									 NoOfMatrices);



			//Rcout <<'\n' <<  i << " retCnts: " << tmpCDr->retCnts.size();
			i++;
			i 			= done;
			for (; i < toWhere; i++){
				rareStruct* CDrAsync;
				CDrAsync = tt[i-done].get();
				// append diversity measures
				divvs.push_back(CDrAsync->div);

				if(NoOfMatrices > 0){
					// append vector to matrix
					uint repI = 0;
					while(repI < CDrAsync->cnts.size()){
						retCnts[repI].push_back(CDrAsync->cnts[repI]);
						repI++;
					}
					// save sample name for naming purposes
					if(CDrAsync->cntsName.size() != 0){
						cntsNames.push_back(CDrAsync->cntsName);
					}
				}
				// skippedNames
				if(CDrAsync->skippedNames.size() > 0){
					skippedNames.push_back(CDrAsync->skippedNames);
				}
				delete CDrAsync;
			}

			// main thread
			divvs.push_back(tmpCDr->div);
			if(NoOfMatrices > 0){
				uint repI = 0;
				while(repI < tmpCDr->cnts.size()){
					retCnts[repI].push_back(tmpCDr->cnts[repI]);
					repI++;
				}
				// save sample name for naming purposes
				if(tmpCDr->cntsName.size() != 0){
					cntsNames.push_back(tmpCDr->cntsName);
				}
			}
			// skippedNames
			if(tmpCDr->skippedNames.size() > 0){
				skippedNames.push_back(tmpCDr->skippedNames);
			}


			delete tmpCDr;



			i++;
			done = i;
		}
		*/
		delete Mo;

		// compute chao2, ACE, ICE and write to file
		computeChao2(chao2, abundInRow);
		computeCE(ICE, abundInRow);
		computeCE(ACE, occuencesInRow);


	}else if(mode == "swap"){
		//rareLowMem(inF, outF, opts->write,  rareDep,  repeats,
		//divvs, retCnts, cntsNames, skippedNames, rowNames, numThr, verbose);
	}


	return 0;
}




