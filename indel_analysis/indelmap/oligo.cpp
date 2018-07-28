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

void loadOligoLookup(std::vector<Oligo*> &output, std::string &exp_oligo_filename) {

	std::string oligo_id, sequence, pam_dir;
	int pam_idx, cut_idx;

	std::ifstream ifs(exp_oligo_filename.c_str(), std::ifstream::in);
	if (!ifs.good()) { std::cout << "Trouble opening expected oligo file " << exp_oligo_filename.c_str(); }
	while (ifs.good()) {
		std::string hline, sline;
		getline(ifs, hline);	//Header line >...
		if (!ifs.good()) break;
		getline(ifs, sequence);	//Sequence line

								//Trim any weird stuff from the end of the sequence
		char last_char = sequence[sequence.length() - 1];
		while (last_char != 'A' && last_char != 'T' && last_char != 'G' && last_char != 'C' && last_char != 'N') {
			sequence = sequence.substr(0, sequence.length() - 1);
			last_char = sequence[sequence.length() - 1];
		}

		//Check for anything crazy remaining
		for (int i = 0; i < sequence.length(); i++) {
			if (sequence[i] != 'A' && sequence[i] != 'T' && sequence[i] != 'G' && sequence[i] != 'C' && sequence[i] != 'N') {
				std::cout << "Invalid character in sequence at position" << i << std::endl;
			}
		}

		if (hline[0] != '>') {
			std::cout << "Invalid fasta input...expecting line starting with >" << std::endl;
			break;
		}
		std::stringstream ss1(hline);
		ss1 >> oligo_id >> pam_idx >> pam_dir;
		oligo_id = oligo_id.substr(1, oligo_id.size() - 1);
		if (pam_dir == "FORWARD") cut_idx = pam_idx - 3;
		else cut_idx = pam_idx + 2;
		if (cut_idx >= sequence.length() - 4) {
			std::cout << oligo_id << " cut site not within sequence (must be before 4 from the end), moving from " << cut_idx << " to " << sequence.length() - 5 << std::endl;
			cut_idx = sequence.length() - 5;
		}
		if (cut_idx < 4) {
			std::cout << oligo_id << " cut site not within sequence (must be at least 4 from the start), moving from " << cut_idx << " to 4" << std::endl;
			cut_idx = 4;
		}
		output.push_back(new Oligo(oligo_id, sequence, cut_idx, pam_dir != "FORWARD"));

	}
	ifs.close();
}

bool subsumedMH(std::vector<Microhomology> &mhs, int i, int j, int mh_len) {
	std::vector<Microhomology>::iterator it = mhs.begin();
	for (; it != mhs.end(); ++it) {
		int delta = it->mh_len - mh_len;
		if (delta <= 0) continue;
		if (it->left_start == (i - delta) && it->right_start == (j - delta)) 
			return true;
	}
	return false;
}

void findAllMicrohomologies(std::vector<Microhomology> &output_mhs, Oligo *oligo, int max_cut_dist) {
	// Collect all
	int N = oligo->seq.length();
	for (int i = 0; i < std::min(N, oligo->cut_idx + max_cut_dist-1); i++) {
		for (int j = std::max(i+1,oligo->cut_idx - max_cut_dist); j < N; j++) {
			int mh_len = 0;
			while ((oligo->seq[i + mh_len] == oligo->seq[j + mh_len]
				     || oligo->seq[i + mh_len] == 'N'
				     || oligo->seq[j + mh_len] == 'N')
				   && (j + mh_len < N)) 
				mh_len++;
			if (mh_len > 1 && !subsumedMH(output_mhs, i, j, mh_len))
				output_mhs.push_back(Microhomology(i,j,mh_len));
		}
	}

}

std::string Microhomology::toString() {
	std::stringstream ss;
	ss << left_start << ":" << right_start << ":" << mh_len;
	return ss.str();
}
