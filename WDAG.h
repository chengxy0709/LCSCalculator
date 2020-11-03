#ifndef WDAG_H
#define WDAG_H

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <list>
#include <limits.h>

#include "point.h"
#include "tool.h"
#include "phash.h"

using namespace std;

class WDAG{

public:	
	WDAG(vector<string>& seqs, string alphabets);
	~WDAG(){}
	
	void run();
	vector<string> MLCS() const {return mlcs;}

private:
	void cal_nums_and_merge_char();
	void cal_weights();
	int  weight(int p, int k);
	void construct_WDAG();
	void determine_levels();
	void GetMLCS(string& LCSRecord, int index);

	vector<string> seqs;
	vector<string> seqms;
	vector<string> mlcs;
	vector< vector< vector<int> > > SucTabs;
	vector< vector<int> > nums;
	vector< vector< vector<int> > > Weights;

	BiHashTable DM;
	vector< list<int> > Edge;
	vector<int> Indegree;
	vector<int> Level;
	vector< list<int> > precursor;

};

#endif //WDAG_H
