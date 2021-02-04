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

#include "gen.h"
#include "readmap.h"

bool getIndelForReadIfMappable(Oligo *oligo, std::string &read_seq, std::string &indel, int max_cut_dist, std::vector<barcode_t> &barcode_lookups) {
	
	std::vector<Oligo*> oligo_lookup; oligo_lookup.push_back(oligo);
	std::vector<barcode_t>::iterator it = barcode_lookups.begin();
	for (; it != barcode_lookups.end(); ++it) {
		buildBarcodeLookup((*it).second, (*it).first, oligo_lookup);
	}
	std::string muts;
	int num_checked, used_barcode;
	std::string oligo_id = getBestOligoIdForRead(read_seq, oligo_lookup, barcode_lookups, indel, muts, true, max_cut_dist, &num_checked, &used_barcode);
	return (oligo_id == oligo->id);
}

void getAllInsertionsOfSize(std::vector<std::string> &inserts, int isize) {
	const char *nts = "ATGC";
	for (int i = 0; i < std::pow(4, isize); i++) {
		std::string ins = "";
		int tmp = i;
		for (int j = 0; j < isize; j++) {
			ins += nts[tmp % 4];
			tmp = tmp / 4;
		}
		inserts.push_back(ins);
	}
}

IndelXX gen_indel_xx;

void generateAllIndels(gen_indel_t &gen_indels, Oligo* oligo, int max_cut_dist, int max_del_size, int max_ins_size, int max_del_to_allow_insert, std::vector<barcode_t> &barcode_lookups, bool check_mappable, rep_reads_t &rep_reads, bool store_reads ) {
	
	if (oligo->reverse) {
		std::cerr << "Warning: Cannot process " << oligo->id << " - REVERSE sequences not supported for indel generation, please take reverse complement and update cut idx." << std::endl;
		return;
	}

	std::string ins_str;
	int left_min = std::max(0, oligo->cut_idx - max_del_size - 1);
	int right_max = std::min(int((oligo->seq).size()),oligo->cut_idx+max_del_size+1);
	for (int left = left_min; left < oligo->cut_idx; left++) {
		for (int right = oligo->cut_idx; right < right_max; right++) {
			int dsize = right - left -1;
			if(dsize > max_del_size) continue;
			int max_ins = max_ins_size * (dsize <= max_del_to_allow_insert);
			for (int isize = 0; isize <= max_ins; isize++) {
				std::vector<std::string> inserts; getAllInsertionsOfSize(inserts, isize);
				std::vector<std::string>::iterator it = inserts.begin();
				for (; it != inserts.end(); ++it) {
					std::string read_seq = oligo->seq.substr(0, left + 1) + *it + oligo->seq.substr(right, (oligo->seq).size() - right);
					if (read_seq.length() < 2) continue;
					std::string indel, muts; bool mappable = true;
					if (check_mappable) {
						mappable = getIndelForReadIfMappable(oligo, read_seq, indel, max_cut_dist, barcode_lookups);
					}else {
						double score = gen_indel_xx.computeMatchScore(*oligo, read_seq, indel, muts, max_cut_dist);
					}
					if (mappable && (indel != "-")) {
						if (gen_indels.find(indel) == gen_indels.end()) {
							gen_indels[indel];
							if (store_reads) rep_reads[indel] = read_seq;
						}
						gen_indels[indel].push_back(loc_t(left, right, *it));
					}
				}
			}
		}
	}
}

char revCmp(char nt) {
	char rc_nt;
	if (nt == 'A') rc_nt = 'T';
	else if (nt == 'T') rc_nt = 'A';
	else if (nt == 'G') rc_nt = 'C';
	else if (nt == 'C') rc_nt = 'G';
	return rc_nt;
}

std::string revCmpSeq(std::string &seq) {
	std::string rc_seq;
	std::string::reverse_iterator rit = seq.rbegin();
	for (; rit != seq.rend(); ++rit) {
		rc_seq += revCmp(*rit);
	}
	return rc_seq;
}

i1_rpt_nts generateI1to3Indels(gen_i1_t &gen_indels, Oligo* oligo, int max_cut_dist) {
	int iloc = oligo->cut_idx + int(oligo->reverse);
	for (int ilen = 1; ilen < 4; ilen++) {
		std::vector<std::string> inserts; getAllInsertionsOfSize(inserts, ilen);
		std::vector<std::string>::iterator it = inserts.begin();
		for (; it != inserts.end(); ++it) {
			std::string read_seq = oligo->seq.substr(0, iloc) + *it + oligo->seq.substr(iloc, (oligo->seq).size());
			std::string indel, muts;
			double score = gen_indel_xx.computeMatchScore(*oligo, read_seq, indel, muts, max_cut_dist);
			if (oligo->reverse) gen_indels[revCmpSeq(*it)] = indel;
			else gen_indels[*it] = indel;
		}
	}
	char left_rpt_nt = oligo->seq[iloc - int(!oligo->reverse)];
	char right_rpt_nt = oligo->seq[iloc - int(oligo->reverse)];
	if (oligo->reverse) {
		left_rpt_nt = revCmp(left_rpt_nt);
		right_rpt_nt = revCmp(right_rpt_nt);
	}
	return i1_rpt_nts(left_rpt_nt, right_rpt_nt);
}

void genIndelsForMicrohomologies(std::vector<std::string> &output_indels, Oligo* oligo, std::vector<Microhomology> &mhs, int max_cut_dist, std::vector<barcode_t> &barcode_lookups) {

	int N = oligo->seq.length();
	std::vector<Microhomology>::iterator mit = mhs.begin();
	for (; mit != mhs.end(); ++mit) {
		std::string exp_read = oligo->seq.substr(0, mit->left_start) + oligo->seq.substr(mit->right_start, N);
		std::string indel; bool mappable = getIndelForReadIfMappable(oligo, exp_read, indel, max_cut_dist, barcode_lookups);
		if(mappable) output_indels.push_back(indel);
		else output_indels.push_back("Unmappable");
	}

}

void writeMhIndelsToFile(std::ofstream &ofs, Oligo *oligo, std::vector<Microhomology> &mhs, std::vector<std::string> &indels) {
	ofs << oligo->id << "\t";
	std::vector<Microhomology>::iterator mit = mhs.begin();
	for (; mit != mhs.end(); ++mit) ofs << mit->toString() << ",";
	ofs << "\t";
	std::vector<std::string>::iterator it = indels.begin();
	for (; it != indels.end(); ++it) ofs << *it << ",";
	ofs << "\n";
}

void writeGeneratedIndelsToFile(std::ofstream &ofs, gen_indel_t &indels, rep_reads_t &rep_reads, int offset) {
	gen_indel_t::iterator git = indels.begin();
	for (; git != indels.end(); ++git) {
		std::vector<loc_t> *locs = &(git->second);
		ofs << git->first << "\t" << locs->size() << "\t";
		std::vector<loc_t>::iterator it = locs->begin();
		ofs << "[" << "(" << it->left+offset << "," << it->right + offset << "," << it->ins_seq << ")"; it++;
		for (; it != locs->end(); ++it)
			ofs << ",(" << it->left + offset << "," << it->right + offset << "," << it->ins_seq << ")";
		ofs << "]";
		if (rep_reads.find(git->first) != rep_reads.end()) ofs << "\t" << rep_reads[git->first];
		ofs << "\n";
	}
}