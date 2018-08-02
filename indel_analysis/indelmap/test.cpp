/*#########################################################################
# Read Mapper for self-target maps
# - Bespoke aligner allowing for crispr edits but otherwise only small (1-2bp)
#   insertions/deletions and mutations.
#########################################################################*/

#include "test.h"
#include "indel.h"
#include "gen.h"

#include <iostream> 
#include <exception>
#include <vector>

Test::Test(){
	description = "None";
}

void Test::run(){
	try{
		excpt_occurred = false;
		passed = false;
		runTest();
	}
	catch(std::exception e){
		excpt_occurred = true;
		std::cout << e.what();
		throw e;
	}
}

void Test::runTest(){
	std::cout << "Running parent runTest" << std::endl;
	//Do nothing: placeholder function
}

void runTests(std::vector<Test*> &tests){

	int num_passed = 0;
	int num_failed = 0;
	int num_run = 0;
	int num_exceptions = 0;

	std::vector<Test*>::iterator it = tests.begin();
	for( ; it != tests.end(); ++it ){
		std::cout << (*it)->description << "...";
		(*it)->run();
		if((*it)->excpt_occurred ){
			num_exceptions++;
			std::cout << "EXCEPTION" << std::endl;
		}else if((*it)->passed ){
			num_passed ++;
			std::cout << "PASS" << std::endl;
		}
		else{ 
			std::cout << "FAIL" << std::endl;
			num_failed++;
		}
		num_run++;
	}

	std::cout << std::endl << "----------------------------------------------" << std::endl;
	std::cout << num_run << " Tests Run" << std::endl;
	std::cout << num_passed << " Passed, ";
	std::cout << num_failed << " Failed, ";
	std::cout << num_exceptions << " Exceptions";
	std::cout << std::endl << "----------------------------------------------" << std::endl;

}

class TestMHIndels : public Test {
public:
	TestMHIndels();
	void runTest();
};

TestMHIndels::TestMHIndels() {
	description = "Test microhomology-related indel generation";
}


class TestGenerateIndels : public Test {
public:
	TestGenerateIndels();
	void runTest();
};

TestGenerateIndels::TestGenerateIndels() {
	description = "Test indel generation";
}


class TestExampleIndels : public Test {
public:
	TestExampleIndels();
	void runTest();
};

TestExampleIndels::TestExampleIndels() {
	description = "Test indel examples";
}

void TestGenerateIndels::runTest() {
	bool pass = true;
	std::string tseq = "AGCTGCGC";
	std::string id = "test";
	Oligo oligo(id, tseq, 4, false);
	int max_cut_dist = 1;
	gen_indel_t gen_indels;
	
	std::vector<barcode_t> barcode_lookups;
	barcode_lookups.push_back(barcode_t(barcode_idx_t(0, 1), barcode_lib_t()));
	barcode_lookups.push_back(barcode_t(barcode_idx_t(-1, 0), barcode_lib_t()));

	generateAllIndels(gen_indels, &oligo, max_cut_dist,5,1,0, barcode_lookups);
	if (gen_indels.size() != 14) {
		std::cout << "Unexpected number of generated indels, expected 36 but found " << gen_indels.size() << std::endl;
		pass = false;
	}
	passed = pass;
}

void TestMHIndels::runTest() {
	bool pass = true;
	std::string tseq = "AATCATGCTCTATGCGAGAC";
	std::string id = "test";
	Oligo oligo(id, tseq, 8, false);
	int max_cut_dist = 1;
	std::vector<Microhomology> mhs;
	findAllMicrohomologies(mhs, &oligo, max_cut_dist);
	if (mhs.size() != 4 || 
		mhs[0].left_start != 1 || mhs[0].right_start != 11 || mhs[0].mh_len != 2 ||
		mhs[1].left_start != 2 || mhs[1].right_start != 8 || mhs[1].mh_len != 2 ||
		mhs[2].left_start != 4 || mhs[2].right_start != 11 || mhs[2].mh_len != 4 ||
		mhs[3].left_start != 7 || mhs[3].right_start != 9 || mhs[3].mh_len != 2) {
		std::cout << "Unexpected microhomolies detected" << std::endl;
		pass = false;
	}

	std::vector<barcode_t> barcode_lookups;
	barcode_lookups.push_back(barcode_t(barcode_idx_t(0, 1), barcode_lib_t()));
	barcode_lookups.push_back(barcode_t(barcode_idx_t(-1, 0), barcode_lib_t()));

	std::vector<std::string> indels;
	genIndelsForMicrohomologies(indels, &oligo, mhs, max_cut_dist, barcode_lookups);
	if (indels.size() != 4 ||
		indels[0] != "D10_L-8C2R5" || indels[1] != "D6_L-7C2R2" ||
		indels[2] != "D7_L-5C4R7" || indels[3] != "D2_L-2C2R3") {
		std::cout << "Unexpected indels detected" << std::endl;
		pass = false;
	}

	passed = pass;
}

class IndelTestCase {
public:
	IndelTestCase(const char *id, const char *tseq, int cut_idx, bool reverse, const char *rseq, const char *exp_indel, const char *exp_mut, int exp_score) :
		id(id), template_seq(tseq), cut_idx(cut_idx), reverse(reverse), read_seq(rseq), exp_indel(exp_indel), exp_muts(exp_mut), exp_score(exp_score) {};
	std::string read_seq, exp_indel, exp_muts;
	Oligo getOligo() { return Oligo(id, template_seq, cut_idx, reverse); };
	int exp_score;
private:
	std::string template_seq, id;
	int cut_idx;
	bool reverse;
};

void TestExampleIndels::runTest() {

	bool pass = true;
	IndelXX indel_xx;

	std::vector<IndelTestCase> test_cases;
	test_cases.push_back(IndelTestCase("Test 22", "CCTCACCTTCCCACGCACGCAGACTCGCAGACGCCCTCTGCTAGAACTGACACGCAGACATTCAGCGGCTCCGACGAG", 38, false, "CCTCACCTTCCCACGCACGCAGACTCGCAGACGCTAGAACTGACACGCAGACATTCAGCGGCTCCGACGAG", "D7_L-7C2R3", "-", 142));
	test_cases.push_back(IndelTestCase("Test 21", "CCTCACCTTCCCACGCACGCAGACTCGCAGACGCCCTCTGCTAGAACTGACACGCAGACATTCAGCGGCTCCGACGAG", 38, false, "CCTCACCTTCCCACGCACGCAGACTCGCAGACGCTGGAACTGACACGCAGACATTCAGCGGCTCCGACGAG", "D7_L-7C2R3", "M4", 139));
	test_cases.push_back(IndelTestCase("Test 20", "AGAAGAGGTTGAAGTGTGCATTGCCACCTCAGTGGTACTTAATGCTGATATGCTCACAATTGCT", 26, true, "AGAAGAGGTTGAAGTGTGCATTGCCACACCAGTGGTACTTAATGCTGATATGCTCACAATTGCT", "-", "M-2M-1", 124));
	test_cases.push_back(IndelTestCase("Test 19", "CAGATGCGATGACCTTTGTGTACTGCTTGTAGCTGGCTGGCCTGGCTGTCATGGCAGTGGGCATCTGGACGCTGGCCCTCAAGAGTGAC", 80, true, "CAGATGCGATGACCTTTGTGTACTGCTTGTAGCTGGCTGGCCTGGCTGTCATGGCAGTGGGCATCTGGACGCTGGCCCTCAAGAGTGA", "-", "D1S-8", 168));
	test_cases.push_back(IndelTestCase("Test 18", "CAGATGCGATGACCTTTGTGTACTGCTTGTAGCTGGCTGGCCTGGCTGTCATGGCAGTGGGCATCTGGACGCTGGCCCTCAAGAGTGAC", 80, true, "AGATGCGATGACCTTTGTGTACTGCTTGTAGCTGGCTGGCCTGGCTGTCATGGCAGTGGGCATCTGGACGCTGGCCCTCAAGAGTGAC", "-", "D1S80", 168));
	test_cases.push_back(IndelTestCase("Test 17", "CAGATGCGATGACCTTTGTGTACTGCTTGTAGCTGGCTGGCCTGGCTGTCATGGCAGTGGGCATCTGGACGCTGGCCCTCAAGAGTGAC", 80, true, "CAGATGCGATGACCTTTGTGTACTGCTTGTAGCTGGCTGGCCTGGCTGTCATGGCAGTGGGCATCTGGACGCTGGCCCTCAAGAGTGACT", "-", "I1S-9", 170));
	test_cases.push_back(IndelTestCase("Test 16", "CAGATGCGATGACCTTTGTGTACTGCTTGTAGCTGGCTGGCCTGGCTGTCATGGCAGTGGGCATCTGGACGCTGGCCCTCAAGAGTGAC", 80, true, "TCAGATGCGATGACCTTTGTGTACTGCTTGTAGCTGGCTGGCCTGGCTGTCATGGCAGTGGGCATCTGGACGCTGGCCCTCAAGAGTGAC", "-", "I1S81", 170));
	test_cases.push_back(IndelTestCase("Test 1", "AAACAGTGCGTCGGCAGGTGCGGATCCTTTCTCCCGACATGTCCGGGATCTGCTACGGGAAGTCAGCCATTTGTTTGTG", 53, false, "AAACAGTGCGTCGGCAGGTGCGGATCCTTTCTCCCGACATGTCCGGGATAGAAGTCAGCCATTTGTTTGTG", "D9_L-5I1R5", "-", 140));
	test_cases.push_back(IndelTestCase("Test 2", "GAATGTTTACCGTGCATTCCCGCTTCATGTATCCTCGCCCGTAGGAGGTATCAAGGATTAGCCACGGAACAATCACGTG", 39, false, "GAATGTTTACCGTGCATTCCACTCATGTATCCTCGCCCGTAGGAGGTATCAAGGATTAGCCCGGAACAATCACGTG", "-", "M-16M-17D2S-18D1S24", 130));
	test_cases.push_back(IndelTestCase("Test 3", "AAACAGTGCGTCGGCAGGTGCGGATCCTTTCTCCCGACATGTCCGGGATCTGCTACGGGAAGTCAGCCATTTGTTTGTG", 53, false, "AAACAGTGCGTCGGCAGGTGCGGATCCTTTCTCCCGACATGATCCTTTCTCCCGACATGTACGGGAAGTCAGCCATTTGTTTGTG", "I18_L-13D12C3R0", "-", 134));
	test_cases.push_back(IndelTestCase("Test 4", "AAATTGGAGTATGGTCACAGTATACTCCCGTAGATTCGGCGGCGGCGACCTGGCGTCTTGGAATAATATGGTACGTGTT", 39, false, "AAATTGGAGTATGGTCACAGTATACTCCCGTAGATTCGGGGACCTGGCGTCTTGGAATAATATGGTACGTGTT", "D7_L-1I1R7", "-", 144));
	test_cases.push_back(IndelTestCase("Test 5", "AAATTGGAGTAGTACAGTACTCCCGTAGATTCGGCGGCGGCGACCTGGCGTCTTGGAATAATATGGTACGTGTT", 33, false, "AAATTGGAGTATGGTCACAGTATACTCCCGTAGATTCGGGGACCTGGCGTCTTGGAATAATATGGTACGTGTT", "D7_L0I1R8", "I2S-17I1S-21I2S-23", 110));
	test_cases.push_back(IndelTestCase("Test 6", "TAAGTTTAATGGATCCAAAAAGTTATACTTAGACCTACAAAAGCACCCATATTACGGGGCACCAGAAGTTGCGGAAGTA", 52, true, "TAAGTTTAATGGATCCAAAAAGTTATACTTAGACCTACAAAAGCACCCATATTACGGGGCACCGAAGATGCGGAAGTA", "-", "D1S-11M-16", 145));
	test_cases.push_back(IndelTestCase("Test 7", "TAAGTTTAATGGATCCAAAAAGTTATACTTAGACCTACAAAAGCACCCATATTACGGGGCACCAGAAGTTGCGGAAGTA", 52, true, "TAAGTTTAATGGATCCAAAAAGTTATACTTAGACCTACAAAAGCACCCATAAGTTGCGGAAGTA", "D15_L-14C1R3", "-", 128));
	test_cases.push_back(IndelTestCase("Test 8", "NAGCNNGCTAGCTNGCTAGCTAGCNAGCTAN", 10, false, "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG", "-", "N-5[A]N-6[T]N-10[T]N3[A]N14[T]N20[G]", 56));
	test_cases.push_back(IndelTestCase("Test 9", "TAAGTTTAATGGATCCAAAAAGNTATNNTTANACCTACAAAAGCACCCATATTACNGGGCACCAGAAGTTGCGGAAGTA", 52, true, "TAAGTTTAATGGATCCAAAAAGTTATACTTAGACCTACAAAAGCACCCATATTACGGGGCACCAGAAGTTGCGGAAGTA", "-", "N21[G]N25[C]N26[A]N30[T]N-3[G]", 153));
	test_cases.push_back(IndelTestCase("Test 10", "ACTAGCTATCAATCGATCGATCGTAGCGAC", 10, false, "ACTAGCTATCGATCGATCGATCGTAGCGAC", "-", "M0", 58));
	test_cases.push_back(IndelTestCase("Test 11", "ACTAGCTATCTATCGATCGATCGTAGCGAC", 10, true, "ACTAGCTATCGATCGATCGATCGTAGCGAC", "-", "M0", 58));
	test_cases.push_back(IndelTestCase("Test 12", "ACTAGCTATCTATCGATCGATCGTAGCGAC", 10, false, "ACTAGCTATCGTCCGATCGATCGTAGCGAC", "-", "M0M1M2", 54));
	test_cases.push_back(IndelTestCase("Test 13", "ACTAGCTATCTATCGATCGATCGTAGCGAC", 10, true, "ACTAGCTATCGTCCGATCGATCGTAGCGAC", "-", "M-2M-1M0", 54));
	test_cases.push_back(IndelTestCase("Test 14", "ACATTTCTGTTACAGTCCCTACCTCATGTGTCAATTCCCTACTGGACAGCACCACCATTTCTAAGGCGGGATGTCTACG", 39, false, "ACATTTCTGTTACAGTCCCTACCTCATGTGTCAATTCCCTACTGGACAGCACCACCATTTCTAAGGCGGGATGTCTACA", "-", "M39", 155));
	test_cases.push_back(IndelTestCase("Test 15", "AAGCCACTTGCGAATAAACCCTGAATTATGGGACAATCCTATGGGACAATTGTGCAGCATCAAAAAGCGGCCAGCTTTT", 39, false, "CAAGGACAATCCTATGGGACAATTGTGCAGCATCAGAAAGCGGCCAGCTTTT", "D28_L-33I1R-4", "M-35M-36M-37M-39M23", 87));
	test_cases.push_back(IndelTestCase("Test 15", "TGGCACTTCTATGCTTCGGTTAGAGAGACAAATCATGATCGCTCATCCTCTACCGTGTTTGCTCTCGGTACCGCGGCAT", 51, true, "TGGCACTTCTATGCTTCGGTTAGAGAGACAAATCATGATACTCATCCGGCAT", "D26_L-23C1R5", "M11D1S12", 93));
	
	std::vector<IndelTestCase>::iterator it = test_cases.begin();
	std::cout << std::endl;
	for (; it != test_cases.end(); ++it) {
		Oligo oligo = it->getOligo();
		std::string indel, muts;
		std::cout << oligo.id << "...";
		int score = indel_xx.computeMatchScore(oligo, it->read_seq, indel, muts, 3);
		if (indel != it->exp_indel) {
			std::cout << oligo.id << " Unexpected Indel: " << indel << " vs " << it->exp_indel << std::endl;
			pass = false;
		}
		if (muts != it->exp_muts) {
			std::cout << oligo.id << " Unexpected Mutations: " << muts << " vs " << it->exp_muts << std::endl;
			pass = false;
		}
		if (score != it->exp_score) {
			std::cout << oligo.id << " Unexpected Score: " << score << " vs " << it->exp_score << std::endl;
			pass = false;
		}
		std::cout << "Done" << std::endl;
	}
	passed = pass;
}

int main(int argc, char *argv[])
{
	std::vector<Test*> tests;
	tests.push_back( new TestMHIndels() );
	tests.push_back( new TestGenerateIndels() );
	tests.push_back( new TestExampleIndels() );
	runTests( tests );
	std::vector<Test*>::iterator it = tests.begin();
	for (; it != tests.end(); ++it) delete (*it);
	return(0);    
}



