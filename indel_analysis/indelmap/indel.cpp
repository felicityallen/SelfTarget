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

void IndelXX::printYStart() {
	std::cout << std::endl;
	for (int j = 0; j < 10; j++) {
		for (int i = 0; i < 10; i++) {
			std::cout << y_left[i][j][0][0] << " ";
		}
		std::cout << std::endl;
	}

}

void IndelXX::printYEnd() {
	std::cout << std::endl;
	for (int j = MAX_R+1; j > MAX_R-10; j--) {
		for (int i = MAX_T+1; i > MAX_T - 10; i--) {
			std::cout << y_right[i][j][0][0] << " ";
		}
		std::cout << std::endl;
	}

}

void IndelXX::initY() {
	y_left[0][0][0][0] = 0; y_left[0][0][0][1] = -1000;
	for (int len = 1; len <= MAX_INDEL; len++) {
		y_left[0][0][len][1] = -1000;
		y_left[0][0][len][0] = -1000;
	}
	y_right[MAX_T+1][MAX_R+1][0][0] = 0; y_right[MAX_T + 1][MAX_R + 1][0][1] = -1000;
	for (int len = 1; len <= MAX_INDEL; len++) {
		y_right[MAX_T + 1][MAX_R + 1][len][1] = -1000;
		y_right[MAX_T + 1][MAX_R + 1][len][0] = -1000;
	}

	for (int i = 1; i < MAX_T + 2; i++) { 
		for (int len = 0; len <= MAX_INDEL; len++) {
			y_left[i][0][len][0] = -GAP_PENALTY;  y_left[i][0][len][1] = -GAP_PENALTY;
		}
	}
	for (int j = 1; j < MAX_R + 2; j++) { 
		for (int len = 0; len <= MAX_INDEL; len++) {
			y_left[0][j][len][0] = -GAP_PENALTY; y_left[0][j][len][1] = -GAP_PENALTY;
		}
	}
	for (int j = 0; j < MAX_R + 1; j++) {
		for (int len = 0; len <= MAX_INDEL; len++) {
			y_right[MAX_T + 1][j][len][0] = -GAP_PENALTY; y_right[MAX_T + 1][j][len][1] = -GAP_PENALTY;
		}
	}
	
	for (int i = 0; i < MAX_T + 1; i++) {
		for (int len = 0; len <= MAX_INDEL; len++) {
			y_right[i][MAX_R+1][len][0] = -GAP_PENALTY;  y_right[i][MAX_R+1][len][1] = -GAP_PENALTY;
		}
	}
}

std::string to_string(int val);

std::string to_string(int val) {
	std::stringstream ss;
	ss << val;
	std::string str = ss.str();
	return str;
}

void IndelXX::populateY( Oligo &oligo, std::string &read_seq, bool left, int max_cut_dist) {

	int T = oligo.seq.length();
	int R = read_seq.length();

	//Set the ranges
	int (*y_ptr)[MAX_T + 2][MAX_R + 2][MAX_INDEL + 1][2] = &y_left;
	int start_i = 1; int end_i = oligo.cut_idx + max_cut_dist + 1; int adj = 1;
	int start_j = 1; int end_j = R + 1; int T_off = 0; int R_off = 0;
	if (!left) {
		T_off = MAX_T - T; R_off = MAX_R - R;
		start_i = MAX_T; end_i = T_off + oligo.cut_idx- max_cut_dist - 1; adj = -1;
		start_j = MAX_R; end_j = R_off;
		y_ptr = &y_right;
	}

	for( int i = start_i; i != end_i; i += adj){
		for( int j = start_j; j != end_j; j += adj){

			//Non-Gap matches/mismatches

			int match = 2*(oligo.seq[i - 1-T_off] == read_seq[j - 1 - R_off]) - (oligo.seq[i - 1-T_off] != read_seq[j - 1 - R_off]);
			if (oligo.seq[i - 1 - T_off] == 'N') match = 1;
			(*y_ptr)[i][j][0][0] = -1000;
			if ((*y_ptr)[i - adj][j - adj][0][0] + match >(*y_ptr)[i][j][0][0])
				(*y_ptr)[i][j][0][0] = (*y_ptr)[i - adj][j - adj][0][0] + match;
			for (int dlen = 1; dlen <= MAX_INDEL; dlen++) {
				if ((*y_ptr)[i - adj][j - adj][dlen][0] + match > (*y_ptr)[i][j][0][0])
					(*y_ptr)[i][j][0][0] = (*y_ptr)[i - adj][j - adj][dlen][0] + match;	//Finish 1-2bp deletion
			}
			for (int ilen = 1; ilen <= MAX_INDEL; ilen++) {
				if ((*y_ptr)[i - adj][j - adj][ilen][1] + match > (*y_ptr)[i][j][0][0])
					(*y_ptr)[i][j][0][0] = (*y_ptr)[i - adj][j - adj][ilen][1] + match;	//Finish 1-2bp insertion
			}

			//Deletions
			(*y_ptr)[i][j][1][0] = (*y_ptr)[i - adj][j][0][0] - GAP_PENALTY;	
			for (int dlen = 2; dlen <= MAX_INDEL; dlen++)
				(*y_ptr)[i][j][dlen][0] = (*y_ptr)[i - adj][j][dlen - 1][0];

			//Insertions
			(*y_ptr)[i][j][1][1] = (*y_ptr)[i][j - adj][0][0] - GAP_PENALTY;
			for (int ilen = 2; ilen <= MAX_INDEL; ilen++)
				(*y_ptr)[i][j][ilen][1] = (*y_ptr)[i][j - adj][ilen - 1][1];

		}
	}
}

void IndelXX::collectMutations(std::string &muts, Oligo &oligo, std::string &read_seq, bool left, int start_i, int start_j) {

	int adj = 1, end_i = 0, end_j = 0, T_off = 0, R_off = 0;
	int T = oligo.seq.length();
	int R = read_seq.length();
	int(*y_ptr)[MAX_T + 2][MAX_R + 2][MAX_INDEL + 1][2] = &y_left;
	if (!left) {
		y_ptr = &y_right; adj = -1; end_i = MAX_T + 1; end_j = MAX_R + 1;
		T_off = MAX_T - T; R_off = MAX_R - R;
		
	}
	int i = start_i, j = start_j, dlen=0, ilen=0;
	while (i != end_i && j != end_j) {
		int score = (*y_ptr)[i][j][dlen][ilen];
		if (dlen == 0 && ilen == 0) {

			int templateN = 0;
			int match = 2*(oligo.seq[i - 1 - T_off] == read_seq[j - 1 - R_off]) - (oligo.seq[i - 1 - T_off] != read_seq[j - 1 - R_off]);
			if (oligo.seq[i - 1 - T_off] == 'N') { match = 1; templateN = 1; }
			int best_prev_score = -1000, prev_i = i, prev_j = j, prev_ilen = ilen, prev_dlen = dlen;
			if ((*y_ptr)[i - adj][j - adj][0][0] > best_prev_score){
				prev_i = i - adj; prev_j = j - adj; prev_ilen = 0; prev_dlen = 0;
				best_prev_score = (*y_ptr)[prev_i][prev_j][0][0];
			}
			for (int len = 1; len <= MAX_INDEL; len++) {
				if ((*y_ptr)[i - adj][j - adj][len][0] > best_prev_score && (*y_ptr)[i - adj][j - adj][len][0] == (score - match)) {
					prev_i = i - adj; prev_j = j - adj; prev_ilen = 0; prev_dlen = len;
					best_prev_score = (*y_ptr)[prev_i][prev_j][len][0];
				}
				if ((*y_ptr)[i - adj][j - adj][len][1] > best_prev_score && (*y_ptr)[i - adj][j - adj][len][1] == (score - match)) {
					prev_i = i - adj; prev_j = j - adj; prev_ilen = len; prev_dlen = 0;
					best_prev_score = (*y_ptr)[prev_i][prev_j][len][1];
				}
			}
			int loc = prev_i - oligo.cut_idx - 1 - T_off;
			int mut_loc = left*(loc + 1) + !left*(loc - 1);
			if (oligo.reverse) {
				loc = oligo.cut_idx + 1 - prev_i + T_off;
				mut_loc = left*(loc - 1) + !left*(loc + 1);
			}
			if (score < best_prev_score) muts += "M" + to_string(mut_loc);
			if (templateN) {
				muts += "N" + to_string(mut_loc) + "[" + read_seq[j - 1 - R_off] + "]";
			}
			if (prev_dlen > 0) muts += "D" + to_string(prev_dlen) + "S" + to_string(loc);
			if (prev_ilen > 0) muts += "I" + to_string(prev_ilen) + "S" + to_string(loc);
			i = prev_i; j = prev_j; dlen = prev_dlen;  ilen = prev_ilen;
		}
		else if (dlen > 0) {
			i -= adj;
			dlen -= 1;
		}
		else if (ilen > 0) {
			j -= adj;
			ilen -= 1;
		}
	}
	if (i != end_i) {
		if (oligo.reverse) muts += 'D' + to_string(abs(i - end_i)) + "S" + to_string(oligo.cut_idx + 1 - i + T_off);
		else muts += 'D' + to_string(abs(i - end_i)) + "S" + to_string(i - oligo.cut_idx - 1 - T_off);
	}
	if (j != end_j) {
		if (oligo.reverse) muts += 'I' + to_string(abs(j - end_j)) + "S" + to_string(oligo.cut_idx + 1 - i + T_off);
		else muts += 'I' + to_string(abs(j - end_j)) + "S" + to_string(i - oligo.cut_idx - 1 - T_off);
	}
}



void IndelXX::populateSummary( int (*y_sum_ptr)[MAX_R+1][2], bool left, int R, int T, int cut_idx, int max_cut_dist){
	
	int start_i = 1, end_i = cut_idx - max_cut_dist, start_j = 1, end_j = R + 1;
	int(*y_ptr)[MAX_T + 2][MAX_R + 2][MAX_INDEL + 1][2] = &y_left;
	int T_off = 0, R_off = 0;
	if (!left) { 
		T_off = MAX_T - T; R_off = MAX_R - R;
		y_ptr = &y_right; start_i = cut_idx + max_cut_dist + 1 + T_off; end_i = MAX_T+1;
		start_j = R_off+1; end_j = MAX_R + 1;
	}
	for (int j = start_j; j < end_j; j++) {
		int score = -1000, best_i = -1;
		for (int i = start_i; i < end_i; i++) {
			if ((*y_ptr)[i][j][0][0] > score) {
				score = (*y_ptr)[i][j][0][0];
				best_i = i;
			}
		}
		(*y_sum_ptr)[j - R_off][0] = score;
		(*y_sum_ptr)[j - R_off][1] = best_i - T_off;		
	}
}

int IndelXX::updateMinMaxIJ(int left_i, int left_j, int right_i, int right_j, int score, int best_score) {
	if (score == best_score) {
		if (left_i < min_left_i)  min_left_i = left_i;
		if (left_j < min_left_j)  min_left_j = left_j;
		if (left_i > max_left_i) {
			max_left_i = left_i;
		}
		if (left_j > max_left_j)  max_left_j = left_j;
		if (right_i < min_right_i)  min_right_i = right_i;
		if (right_j < min_right_j)  min_right_j = right_j;
		if (right_i > max_right_i)  max_right_i = right_i;
		if (right_j > max_right_j)  max_right_j = right_j;
	}
	else if (score > best_score) {
		min_left_i = left_i; max_left_i = left_i;
		min_left_j = left_j; max_left_j = left_j;
		min_right_i = right_i; max_right_i = right_i;
		min_right_j = right_j; max_right_j = right_j;
		best_score = score;
	}
	return best_score;
}

//Dynamic program for computing a string match that is specific to CRISPR mutations
//(limits large insertions/deletions to those spanning the cut site)
int IndelXX::computeMatchScore( Oligo &oligo, std::string &read_seq, std::string &output_indel_str, std::string &muts, int max_cut_dist){
		
	int T = oligo.seq.length(); int T_off = MAX_T - T;
	int R = read_seq.length(); int R_off = MAX_R - R;
	if (T_off < 0) {
		std::cout << "Template oligo length " << T << " greater than maximum allowed " << MAX_T << std::endl;
		return -1;
	}
	if (R_off < 0) {
		std::cout << "Read length " << R << " greater than maximum allowed " << MAX_R << std::endl;
		return -1;
	}

	//Populate the left and right arrays
	populateY(oligo, read_seq, true, max_cut_dist);
	populateY(oligo, read_seq, false, max_cut_dist);

	//Find the maximum score combining the left and right sides and allowing max_cut_dist bp ambiguity in the template cut site
	int best_score = -1;
	
	//Summmarise the regions away from the cut site
	int y_left_summary[MAX_R+1][2]; populateSummary(&y_left_summary, true, R, T, oligo.cut_idx, max_cut_dist);
	int y_right_summary[MAX_R+1][2]; populateSummary(&y_right_summary, false, R, T, oligo.cut_idx, max_cut_dist);

	//Check for maximum score away from the cut site
	for (int left_j = 1; left_j <= R; left_j++) {
		for (int right_j = left_j + 1; right_j <= R; right_j++) {
			int score = y_left_summary[left_j][0] + y_right_summary[right_j][0];
			int left_i = y_left_summary[left_j][1];
			int right_i = y_right_summary[right_j][1];
			best_score = updateMinMaxIJ(left_i, left_j, right_i, right_j, score, best_score);
		}
	}

	//Check for better scores, with left_i left of the ambiguous region, right_i in the ambiguous region
	for (int right_i = oligo.cut_idx - max_cut_dist; right_i <= oligo.cut_idx + max_cut_dist; right_i++) {
		for (int left_j = 1; left_j <= R; left_j++) {
			int left_i = y_left_summary[left_j][1];
			for (int right_j = left_j+1; right_j <= R; right_j++) {
				int score = y_left_summary[left_j][0] + y_right[right_i+T_off][right_j+R_off][0][0];
				best_score = updateMinMaxIJ(left_i, left_j, right_i, right_j, score, best_score);
			}
		}
	}

	//Check for better scores, with left_i in the ambiguous region, right_i right of the ambiguous region
	for (int right_j = 1; right_j <= R; right_j++) {
		int right_i = y_right_summary[right_j][1];
		for (int left_i = oligo.cut_idx - max_cut_dist; left_i <= oligo.cut_idx + max_cut_dist; left_i++) {
			for (int left_j = 1; left_j < right_j; left_j++) {
				int score = y_right_summary[right_j][0] + y_left[left_i][left_j][0][0];
				best_score = updateMinMaxIJ(left_i, left_j, right_i, right_j, score, best_score);
			}
		}
	}

	//Check for better scores within the ambiguous region around the cut site
	for (int left_i = oligo.cut_idx - max_cut_dist; left_i <= oligo.cut_idx + max_cut_dist; left_i++) {
		for (int right_i = left_i + 1; right_i <= oligo.cut_idx + max_cut_dist; right_i++) {
			for (int left_j = 1; left_j <= R; left_j++) {
				for (int right_j = left_j + 1; right_j <= R; right_j++) {
					int score = y_left[left_i][left_j][0][0] + y_right[right_i+T_off][right_j+R_off][0][0];
					best_score = updateMinMaxIJ(left_i, left_j, right_i, right_j, score, best_score);
				}
			}
		}
	}

	//Don't bother processing the indel if the score is too low
	if (best_score < 0.5*read_seq.length()) return best_score;

	//Interpret what happened
	int del_size = min_right_i - min_left_i - 1;
	int ins_size = min_right_j - min_left_j - 1;
	int central_dsize = max_left_i - min_left_i;
	int central_isize = max_left_j - min_left_j;

	int left = min_left_i - oligo.cut_idx - 1;
	int right = max_right_i - oligo.cut_idx - 1;
	if (oligo.reverse) {
		left = oligo.cut_idx + 1 - max_right_i;
		right = oligo.cut_idx + 1 - min_left_i;
	}
	if (del_size == ins_size && del_size < 4) {
		for (int i = 1; i <= del_size; i++)
			muts += "M" + to_string(left + i);
		output_indel_str = "-";
	}
	else if (del_size >= ins_size && del_size > 0) {
		output_indel_str = "D" + to_string(del_size) + "_L" + to_string(left);
		if (ins_size > 0) 	output_indel_str += "I" + to_string(ins_size);
		if (central_dsize > 0) output_indel_str += "C" + to_string(central_dsize);
		//if (central_isize > 0) output_indel_str += "X" + to_string(central_isize);
		output_indel_str += "R" + to_string(right);
	}
	else if (ins_size > del_size && ins_size > 0) {
		output_indel_str = "I" + to_string(ins_size) + "_L" + to_string(left);
		if (del_size > 0) 	output_indel_str += "D" + to_string(del_size);
		if (central_dsize > 0) output_indel_str += "C" + to_string(central_dsize);
		//if (central_isize > 0) output_indel_str += "X" + to_string(central_isize);
		output_indel_str += "R" + to_string(right);
	}
	else output_indel_str = "-";

	collectMutations(muts, oligo, read_seq, true, min_left_i, min_left_j);
	collectMutations(muts, oligo, read_seq, false, max_right_i + T_off, max_right_j +R_off);
	if (muts.length() == 0) muts = "-";
	
	return best_score;

}
