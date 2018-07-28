/*#########################################################################
# Indel Generater for self-target maps
# - Generates all possible indels for a given oligo target
#########################################################################*/

int main(int argc, char *argv[]);

static const int MAX_DEL_SIZE = 30;
static const int MAX_INS_SIZE = 2;
static const int MAX_DEL_TO_ALLOW_INS = 0;

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <ctime>

#include "oligo.h"
#include "gen.h"
#include "indel.h"
#include "version.h"

void writeGeneratedIndelsToFile(std::ofstream &ofs, gen_indel_t &indels) {
	gen_indel_t::iterator git = indels.begin();
	for (; git != indels.end(); ++git ){
		std::vector<loc_t> *locs = &(git->second);
		ofs << git->first << "\t" << locs->size() << "\t";
		std::vector<loc_t>::iterator it = locs->begin();
		ofs << "[" << "(" << it->left << "," << it->right << "," << it->ins_seq << ")"; it++;
		for (; it != locs->end(); ++it) 
			ofs << ",(" << it->left << "," << it->right << "," << it->ins_seq << ")";
		ofs << "]\n";
	}
}

int main(int argc, char *argv[])
{
	if (argc != 3 && argc != 4)
	{
		std::cout << std::endl << "Usage:" << std::endl;
		std::cout << "indelgen.exe <exp_oligo_file> <output_prefix> <(opt) max_cut_dist (default 4)>" << std::endl << std::endl;
		exit(1);
	}
    std::string exp_oligo_file = argv[1];
	std::string output_prefix = argv[2];
	int max_cut_dist = 4;
	if(argc >= 4) max_cut_dist = atoi(argv[3]);

	//Read in the expected oligos
	std::vector<Oligo*> oligo_lookup;
	loadOligoLookup(oligo_lookup, exp_oligo_file);
	
	//Load the barcode specifications (used to check whether generated indels are detectable)
	std::vector<barcode_t> barcode_lookups; loadDefaultBarcodes(barcode_lookups);

	//Generate possible indels for each oligo and write to file
	std::vector<Oligo*>::iterator it = oligo_lookup.begin();
	for (; it != oligo_lookup.end(); ++it) {
		gen_indel_t indels;
		generateAllIndels(indels, *it, max_cut_dist, MAX_DEL_SIZE, MAX_INS_SIZE, MAX_DEL_TO_ALLOW_INS, barcode_lookups);
		std::string oligo_output_filename = output_prefix + (*it)->id + "_genindels.txt";
		std::ofstream ofs(oligo_output_filename.c_str(), std::fstream::out);
		ofs << "@@@Git Commit: " << GIT_COMMIT_HASH << std::endl;
		writeGeneratedIndelsToFile(ofs, indels);
		ofs.close();
	}
	
	//Free the oligo memory
	std::vector<Oligo*>::iterator ito = oligo_lookup.begin();
	for (; ito != oligo_lookup.end(); ++ito) delete *ito;

	return(0);    
}
