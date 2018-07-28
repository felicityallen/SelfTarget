/*#########################################################################
# Read Mapper for self-target maps
# - Bespoke aligner allowing for crispr edits but otherwise only small (1-2bp) 
#   insertions/deletions and mutations.
#########################################################################*/

#ifndef __OLIGO_H__
#define __OLIGO_H__

#include <vector>
#include <string>
#include <sstream>

class Oligo {
public:
	Oligo(std::string &an_id, std::string &a_seq, int a_cut_idx, bool a_reverse)
	{
		id = an_id; seq = a_seq; cut_idx = a_cut_idx; reverse = a_reverse;
	};
	std::string id;
	std::string seq;
	int cut_idx;
	bool reverse;
};

class Microhomology {
public:
	Microhomology(int a_left_start, int a_right_start, int a_mh_len) :
		left_start(a_left_start), right_start(a_right_start), mh_len(a_mh_len) {};
	int left_start;
	int right_start;
	int mh_len;
	std::string toString();
};

void loadOligoLookup(std::vector<Oligo*> &output, std::string &exp_oligo_filename);
void findAllMicrohomologies(std::vector<Microhomology> &output_mhs, Oligo *oligo, int max_cut_dist);

#endif // __OLIGO_H__