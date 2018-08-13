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
#include <assert.h>

#include "oligo.h"
#include "gen.h"
#include "indel.h"
#include "version.h"

static const int MAX_FLANK = 50;
static const int MIN_FLANK = 10;

int main(int argc, char *argv[])
{
	if (argc != 4)
	{
		std::cout << std::endl << "Usage:" << std::endl;
		std::cout << "indelgentarget.exe <seq> <pam_idx> <output_file>" << std::endl << std::endl;
		exit(1);
	}

	//Single input sequence and cut idx
	int max_cut_dist = 4;
	std::string seq = argv[1];
	int pam_idx = atoi(argv[2]);
	std::string output_file = argv[3];

	//Trim seq if too long
	int left_trim = 0;
	if (seq.size() - pam_idx > MAX_FLANK) {
		seq = seq.substr(0, pam_idx + 50);
	}
	if (pam_idx > MAX_FLANK) {
		seq = seq.substr(pam_idx - MAX_FLANK);
		left_trim = pam_idx - MAX_FLANK;
		pam_idx = MAX_FLANK;
	}

	//Check sequence characteristics
	assert(pam_idx > MIN_FLANK);
	assert(seq.size() - pam_idx > MIN_FLANK);
	assert(seq.substr(pam_idx + 1, 2) == "GG");

	gen_indel_t indels; rep_reads_t rep_reads;  std::vector<barcode_t> barcode_lookups;
	std::string id_str("No Id");
	Oligo oligo(id_str, seq, pam_idx - 3, false);
	generateAllIndels(indels, &oligo, max_cut_dist, MAX_DEL_SIZE, MAX_INS_SIZE, MAX_DEL_TO_ALLOW_INS, barcode_lookups, false, rep_reads, true);
	std::ofstream ofs(output_file.c_str(), std::fstream::out);
	ofs << "@@@Git Commit: " << GIT_COMMIT_HASH << std::endl;
	writeGeneratedIndelsToFile(ofs, indels, rep_reads, left_trim);
	ofs.close();

	return(0);
}
