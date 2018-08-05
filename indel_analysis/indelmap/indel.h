/*#########################################################################
# Read Mapper for self-target maps
# - Bespoke aligner allowing for crispr edits but otherwise only small (1-2bp)
#   insertions/deletions and mutations.
#########################################################################*/

#ifndef __INDEL_H__
#define __INDEL_H__

#include <map>
#include "oligo.h"

static const int MAX_T = 300;
static const int MAX_R = 300;
static const int MAX_INDEL = 2;
static const int GAP_PENALTY = 8;

class IndelXX {
public:
	IndelXX() { initY(); };
	int computeMatchScore(Oligo &oligo, std::string &read_seq, std::string &output_indel_str, std::string &muts, int max_cut_dist);

private:
	void initY();
	void printYStart();
	void printYEnd();
	void populateY(Oligo &oligo, std::string &read_seq, bool left, int max_cut_dist);
	void collectMutations(std::string &muts, Oligo &oligo, std::string &read_seq, bool left, int start_i, int start_j);
	void populateSummary(int(*y_sum_ptr)[MAX_R + 1][2], bool left, int R, int T, int cut_idx, int max_cut_dist);
	int updateMinMaxIJ(int left_i, int left_j, int right_i, int right_j, int score, int best_score);

	int y_left[MAX_T + 2][MAX_R + 2][MAX_INDEL + 1][2];
	int y_right[MAX_T + 2][MAX_R + 2][MAX_INDEL + 1][2];

	int min_left_i, max_left_i, min_right_i, max_right_i;
	int min_left_j, max_left_j, min_right_j, max_right_j;
};

#endif // __INDEL_H__