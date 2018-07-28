/*#########################################################################
# Indel Generater for self-target maps
# - Generates all possible indels for a given oligo target
#########################################################################*/

int main(int argc, char *argv[]);

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
#include "indel.h"
#include "gen.h"
#include "version.h"

int main(int argc, char *argv[])
{
	if (argc != 3 && argc != 4)
	{
		std::cout << std::endl << "Usage:" << std::endl;
		std::cout << "indelmh.exe <exp_oligo_file> <output_file> <(opt) max_cut_dist (default 4)>" << std::endl << std::endl;
		exit(1);
	}
    std::string exp_oligo_file = argv[1];
	std::string output_file = argv[2];
	int max_cut_dist = 4;
	if(argc >= 4) max_cut_dist = atoi(argv[3]);

	//Read in the expected oligos
	std::vector<Oligo*> oligo_lookup;
	loadOligoLookup(oligo_lookup, exp_oligo_file);
	
	//Load the barcode specifications (used to check whether generated indels are detectable)
	std::vector<barcode_t> barcode_lookups; loadDefaultBarcodes(barcode_lookups);

	//Check for microhomologies, and if found, generate matched indels, write to file
	std::ofstream ofs(output_file.c_str(), std::fstream::out);
	ofs << "@@@Git Commit: " << GIT_COMMIT_HASH << std::endl;
	
	std::vector<Oligo*>::iterator it = oligo_lookup.begin();
	for (; it != oligo_lookup.end(); ++it) {
		
		std::vector<Microhomology> mhs;
		findAllMicrohomologies(mhs, *it, max_cut_dist );
		
		std::vector<std::string> indels;
		genIndelsForMicrohomologies(indels, *it, mhs, max_cut_dist, barcode_lookups);

		writeMhIndelsToFile(ofs, *it, mhs, indels);
		
	}
	ofs.close();

	//Free the oligo memory
	std::vector<Oligo*>::iterator ito = oligo_lookup.begin();
	for (; ito != oligo_lookup.end(); ++ito) delete *ito;

	return(0);    
}
