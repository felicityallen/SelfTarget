/*#########################################################################
# Read Mapper for self-target maps
# - Bespoke aligner allowing for crispr edits but otherwise only small (1-2bp) 
#   insertions/deletions and mutations.
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

#include "indel.h"
#include "oligo.h"
#include "readmap.h"
#include "version.h"

int main(int argc, char *argv[])
{
	if (argc != 5 && argc != 6 && argc != 7)
	{
		std::cout << std::endl << "Usage:" << std::endl;
		std::cout << "indelmap.exe <fastq_file> <exp_oligo_file> <output_file> <barcode_check_only (1 or 0)> <(opt) max_cut_dist (default 4)> <(opt)barcode_idx_file>" << std::endl << std::endl;
		exit(1);
	}
    std::string fastq_file = argv[1];
    std::string exp_oligo_file = argv[2];
	std::string output_filename = argv[3];
	int barcode_check_only = atoi(argv[4]);
	int max_cut_dist = 4;
	if(argc >= 6) max_cut_dist = atoi(argv[5]);
	std::vector<barcode_t> barcode_lookups;
	if (argc >= 7) {
		std::string barcode_file = argv[6];
		loadBarcodeIdxs(barcode_lookups, barcode_file);
	}else loadDefaultBarcodes(barcode_lookups);

	//Read in the expected oligos and create barcode hashes
	std::vector<Oligo*> oligo_lookup;
	loadOligoLookup(oligo_lookup, exp_oligo_file);

	//Load the barcode lookups
	std::vector<barcode_t>::iterator it = barcode_lookups.begin();
	for (; it != barcode_lookups.end(); ++it) {
		buildBarcodeLookup((*it).second, (*it).first, oligo_lookup);
	}
	
	//Set up the output
	std::ofstream ofs(output_filename.c_str(), std::fstream::out );
	ofs << "@@@Git Commit: " << GIT_COMMIT_HASH << std::endl;

	//Map the fasta or fastq reads
	bool is_fastq = (fastq_file[fastq_file.length() - 1] == 'q');
	std::ifstream ifs(fastq_file.c_str(), std::ifstream::in);
	if (!ifs.good()) { 
		std::cout << "Trouble opening fastq file " << fastq_file.c_str(); 
	}
	while (ifs.good()) {
		std::string hline, sline, nullline;
		getline(ifs, hline); 	//read id
		if (!ifs.good()) break;
		getline(ifs, sline); 	//Sequence
		if (is_fastq) {	//else Fasta
			if (!ifs.good()) break;
			getline(ifs, nullline); //Ignore 
			if (!ifs.good()) break;
			getline(ifs, nullline);	//Ignore
		}
		if (hline[0] != '@' && hline[0] != '>') {
			std::cout << "Invalid read input...expecting line starting with > or @" << std::endl;
			break;
		}
		clock_t begin_time = clock();

		int num_checked = 0, used_barcode = 0;
		std::string indel, muts;
		std::string oligo_id = getBestOligoIdForRead(sline, oligo_lookup, barcode_lookups, indel, muts, barcode_check_only, max_cut_dist, &num_checked, &used_barcode);
		clock_t end_time = clock();
		double elapsed_secs = double(end_time - begin_time) / CLOCKS_PER_SEC;
		ofs << hline.substr(1,hline.length()) << '\t' << oligo_id << '\t' << indel << '\t' << muts << '\t' << elapsed_secs << "\t" << num_checked << "\t" << used_barcode << std::endl;

	}
	ofs.close();

	//Free the oligo memory
	std::vector<Oligo*>::iterator ito = oligo_lookup.begin();
	for (; ito != oligo_lookup.end(); ++ito) delete *ito;

	return(0);    
}
