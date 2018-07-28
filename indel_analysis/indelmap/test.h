/*#########################################################################
# Read Mapper for self-target maps
# - Bespoke aligner allowing for crispr edits but otherwise only small (1-2bp)
#   insertions/deletions and mutations.
#########################################################################*/

#ifndef __TEST_H__
#define __TEST_H__

#include <string>

class Test{

public:
	Test();
	void run();
	virtual void runTest();
	bool passed;
	bool excpt_occurred;
	std::string description;
	virtual ~Test(){}
};

#endif // __TEST_H__