/*#########################################################################
# Indel Generater for self-target maps
# - Generates all possible indels for a given oligo target
#########################################################################*/

int main(int argc, char *argv[]);

static const int MAX_DEL_SIZE = 0;
static const int MAX_INS_SIZE = 1;

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

char nts[4] = { 'A','T','G','C' };

void writeGeneratedIndelsToFile(std::ofstream &ofs, gen_i1_t &indels, std::string &oligo_id, i1_rpt_nts &rpt_nts) {
	ofs << oligo_id;
	for (int ilen = 1; ilen < 4; ilen++) {
		std::vector<std::string> inserts; getAllInsertionsOfSize(inserts, ilen);
		std::vector<std::string>::iterator it = inserts.begin();
		for (; it != inserts.end(); ++it) {
			ofs << "\t" << indels[*it];
		}
	}
	ofs << "\t" << std::get<0>(rpt_nts) << "\t" << std::get<1>(rpt_nts) << "\n";
}

int main(int argc, char *argv[])
{
	if (argc != 2 && argc != 3)
	{
		std::cout << std::endl << "Usage:" << std::endl;
		std::cout << "indelgen_i1.exe <exp_oligo_file> <(opt) max_cut_dist (default 4)>" << std::endl << std::endl;
		exit(1);
	}
	std::string exp_oligo_file = argv[1];
	int max_cut_dist = 4;
	if (argc >= 3) max_cut_dist = atoi(argv[2]);

	//Read in the expected oligos
	std::vector<Oligo*> oligo_lookup;
	loadOligoLookup(oligo_lookup, exp_oligo_file);

	//Load the barcode specifications (used to check whether generated indels are detectable)
	std::vector<barcode_t> barcode_lookups; loadDefaultBarcodes(barcode_lookups);

	//Generate possible indels for each oligo and write to file
	std::string oligo_output_filename = exp_oligo_file.substr(0,exp_oligo_file.length()-6) + "_gen_i1_indels.txt";
	std::ofstream ofs(oligo_output_filename.c_str(), std::fstream::out);
	ofs << "@@@Git Commit: " << GIT_COMMIT_HASH << std::endl;
	ofs << "Oligo Id";
	for (int ilen = 1; ilen < 4; ilen++) {
		std::vector<std::string> inserts; getAllInsertionsOfSize(inserts, ilen);
		std::vector<std::string>::iterator iit = inserts.begin();
		for (; iit != inserts.end(); ++iit)
			ofs << "\tI1_" << *iit;
	}
	ofs << "\t" << "Repeat Nucleotide Left\tRepeat Nucleotide Right\n";

	std::vector<Oligo*>::iterator it = oligo_lookup.begin();
	for (; it != oligo_lookup.end(); ++it) {
		gen_i1_t indels;
		i1_rpt_nts rpt_nts = generateI1to3Indels(indels, *it, max_cut_dist);
		writeGeneratedIndelsToFile(ofs, indels, (*it)->id, rpt_nts);
	}
	ofs.close();
	//Free the oligo memory
	std::vector<Oligo*>::iterator ito = oligo_lookup.begin();
	for (; ito != oligo_lookup.end(); ++ito) delete *ito;

	return(0);
}
