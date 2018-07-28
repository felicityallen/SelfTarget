/*#########################################################################
# Read Mapper for self-target maps
# - Bespoke aligner allowing for crispr edits but otherwise only small (1-2bp) 
#   insertions/deletions and mutations.
#########################################################################*/

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
#include "readmap.h"

IndelXX readmap_indel_xx;

void loadBarcodeIdxs(std::vector<barcode_t> &barcode_lookups, std::string &barcode_filename) {
	std::ifstream ifs(barcode_filename.c_str(), std::ifstream::in);
	if (!ifs.good()) { std::cout << "Trouble opening expected barcode file " << barcode_filename.c_str(); }
	while (ifs.good()) {
		std::string line;
		int idx1 = -1;
		int idx2 = -2;
		getline(ifs, line);
		std::stringstream ss1(line);
		ss1 >> idx1 >> idx2;
		barcode_lookups.push_back(barcode_t(barcode_idx_t(idx1, idx2), barcode_lib_t()));
	}
	ifs.close();
}


std::string pickBestOligoIdFromList(std::string &read_seq, std::vector<Oligo*> &oligo_lookup, std::string &output_indel, std::string &mutations, int max_cut_dist) {

	std::string best_oligo_id, best_oligo_seq, indel;
	int best_score = read_seq.length()*0.8;	//Needs to map at least 80% of the read sequence

	std::vector<Oligo*>::iterator it = oligo_lookup.begin();
	for (; it != oligo_lookup.end(); ++it) {
		std::string indel, muts;
		double score = readmap_indel_xx.computeMatchScore(**it, read_seq, indel, muts, max_cut_dist);
		if (score == best_score && best_oligo_id.length() != 0) {
			best_oligo_id += ("," + (*it)->id);
			output_indel += ("," + indel);
			mutations += ("," + muts);
		}
		if (score > best_score) {
			best_score = score;
			best_oligo_id = (*it)->id;
			output_indel = indel;
			mutations = muts;
		}
	}
	return best_oligo_id;
}

void insertOligos(std::vector<Oligo*> &input_oligos, std::vector<Oligo*> &input_list) {
	std::vector<Oligo*>::iterator it = input_oligos.begin();
	for (; it != input_oligos.end(); ++it) {
		if (std::find(input_list.begin(), input_list.end(), *it) == input_list.end())
			input_list.push_back(*it);
	}
}

std::string getBarcode(barcode_idx_t bc_idx, std::string &seq) {
	int idx1 = bc_idx.first;
	if (idx1 < 0) idx1 = seq.length() + idx1;
	int idx2 = bc_idx.second;
	if (idx2 <= 0) idx2 = seq.length() + idx2;
	if (idx1 >= seq.length() || idx2 > seq.length()) {
		throw BarcodeNotComputableException();
	}
	return seq.substr(idx1, idx2 - idx1);
}

std::string getBestOligoIdForRead(std::string &read_seq, std::vector<Oligo*> &oligo_lookup, std::vector<barcode_t> &bc_lookups, std::string &output_indel, std::string &mutations, bool barcode_check_only, int max_cut_dist, int *num_checked, int *used_barcode) {
	std::string best_oligo_id;

	//Try to find a match amongst oligos with close barcodes
	std::vector<Oligo*> oligo_list;
	std::vector<barcode_t>::iterator it = bc_lookups.begin();
	for (; it != bc_lookups.end(); ++it) {
		barcode_idx_t idxs = (*it).first;
		barcode_lib_t *bc_lib = &(*it).second;
		try {
			std::string bc = getBarcode(idxs, read_seq);
			if (bc_lib->find(bc) != bc_lib->end())
				insertOligos((*bc_lib)[bc], oligo_list);
		}
		catch (BarcodeNotComputableException) {}
	}
	best_oligo_id = pickBestOligoIdFromList(read_seq, oligo_list, output_indel, mutations, max_cut_dist);
	*num_checked = oligo_list.size();
	*used_barcode = 1;

	//If we found one, return it, else either skip or process full list (based on barcode_check_only)
	if (best_oligo_id.length() > 0) return best_oligo_id;
	else {
		if (barcode_check_only)
			return std::string("None");
		else {
			*used_barcode = 0;
			*num_checked = oligo_lookup.size();
			return pickBestOligoIdFromList(read_seq, oligo_lookup, output_indel, mutations, max_cut_dist);
		}
	}
}

void addOligoToBarcodeLookup(barcode_lib_t &bc_lookup, std::string &bc, Oligo* oligo_ptr) {
	if (bc_lookup.find(bc) == bc_lookup.end())
		bc_lookup[bc];
	if (std::find(bc_lookup[bc].begin(), bc_lookup[bc].end(), oligo_ptr) == bc_lookup[bc].end()) {
		bc_lookup[bc].push_back(oligo_ptr);
	}
}

void addAllBarcodesWithOneOrLessMutations(barcode_lib_t &bc_lookup, std::string &bc, Oligo* oligo_ptr) {
	char nucs[4] = { 'A','T','G','C' };
	for (int i = 0; i < bc.length(); i++) {
		char orig_iN = bc[i];
		for (int n1 = 0; n1 < 4; n1++) {
			bc[i] = nucs[n1]; addOligoToBarcodeLookup(bc_lookup, bc, oligo_ptr);
		}
		bc[i] = orig_iN;
	}
}

void buildBarcodeLookup(barcode_lib_t &bc_lookup, barcode_idx_t idxs, std::vector<Oligo*> &oligo_lookup) {
	std::vector<Oligo*>::iterator it = oligo_lookup.begin();
	for (; it != oligo_lookup.end(); ++it) {
		std::string bc = getBarcode(idxs, (*it)->seq);
		addAllBarcodesWithOneOrLessMutations(bc_lookup, bc, *it);
	}
}

void loadDefaultBarcodes(std::vector<barcode_t> &barcode_lookups) {
	barcode_lookups.push_back(barcode_t(barcode_idx_t(0, 10), barcode_lib_t()));
	barcode_lookups.push_back(barcode_t(barcode_idx_t(5, 15), barcode_lib_t()));
	barcode_lookups.push_back(barcode_t(barcode_idx_t(-15, -5), barcode_lib_t()));
	barcode_lookups.push_back(barcode_t(barcode_idx_t(-10, 0), barcode_lib_t()));
}

