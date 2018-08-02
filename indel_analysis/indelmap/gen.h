/*#########################################################################
# Read Mapper for self-target maps
# - Bespoke aligner allowing for crispr edits but otherwise only small (1-2bp)
#   insertions/deletions and mutations.
#########################################################################*/

#ifndef __GEN_H__
#define __GEN_H__
#include <cmath>
#include <map>
#include "oligo.h"
#include "indel.h"
#include "readmap.h"

typedef struct gen_indel_loc { 
	int left; int right; std::string ins_seq; 
	gen_indel_loc(int a_left, int a_right, std::string a_seq) { left = a_left; right = a_right; ins_seq = a_seq; };
} loc_t;
typedef std::tuple<char, char> i1_rpt_nts;
typedef std::map<std::string, std::vector<loc_t>> gen_indel_t;
typedef std::map<std::string, std::string> gen_i1_t;
void generateAllIndels(gen_indel_t &gen_indels, Oligo* oligo, int max_cut_dist, int max_del_size, int max_ins_size, int max_del_to_allow_insert, std::vector<barcode_t> &barcode_lookups);
void genIndelsForMicrohomologies(std::vector<std::string> &output_indels, Oligo* oligo, std::vector<Microhomology> &mhs, int max_cut_dist, std::vector<barcode_t> &barcode_lookups);
void writeMhIndelsToFile(std::ofstream &ofs, Oligo *oligo, std::vector<Microhomology> &mhs, std::vector<std::string> &indels);
i1_rpt_nts generateI1to3Indels(gen_i1_t &gen_indels, Oligo* oligo, int max_cut_dist);
void getAllInsertionsOfSize(std::vector<std::string> &inserts, int isize);

#endif // __GEN_H__
