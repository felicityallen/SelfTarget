/*#########################################################################
# Read Mapper for self-target maps
# - Bespoke aligner allowing for crispr edits but otherwise only small (1-2bp) 
#   insertions/deletions and mutations.
#########################################################################*/

#ifndef __READMAP_H__
#define __READMAP_H__

#include <vector>
#include <string>
#include <sstream>

#include "oligo.h"

struct BarcodeNotComputableException : public std::exception {
	const char * what() const throw () { return "Could not compute barcode for sequence"; }
};

typedef std::map<std::string, std::vector<Oligo*> > barcode_lib_t;
typedef std::pair<int, int> barcode_idx_t;
typedef std::pair<barcode_idx_t, barcode_lib_t> barcode_t;

void loadBarcodeIdxs(std::vector<barcode_t> &barcode_lookups, std::string &barcode_filename);
void buildBarcodeLookup(barcode_lib_t &bc_lookup, barcode_idx_t idxs, std::vector<Oligo*> &oligo_lookup);
void loadDefaultBarcodes(std::vector<barcode_t> &barcode_lookups);
std::string getBestOligoIdForRead(std::string &read_seq, std::vector<Oligo*> &oligo_lookup, std::vector<barcode_t> &bc_lookups, std::string &output_indel, std::string &mutations, bool barcode_check_only, int max_cut_dist, int *num_checked, int *used_barcode);

#endif // __READMAP_H__