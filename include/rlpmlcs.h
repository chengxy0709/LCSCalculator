#ifndef RLP_MLCS_H
#define RLP_MLCS_H

#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <list>
#include <algorithm>

#include "point.h"
#include "phash.h"
#include "tool.h"

#define ALG_NO_REQUIRED

#define RLPMLCSSYM "RLPMLCS"

class RLP_MLCS{

public:	
	RLP_MLCS(vector<string>& seqs, string& alphabets);
	~RLP_MLCS();
	
	void run();
	vector<string> MLCS() const {return mlcs;}

private:
    int cal_opt_point();
	int construct_ICSG(Point<CordType> *p0);
	int ForwardTopSort(int maxIndex, int startpoint, int curlevel);
	void BackwardTopSort(int maxlevel, int pinf_index);
    void construct_SubICSG(int startpoint, int curlevel);
	void GetMLCS(string& LCSRecord, int index);

	vector<string> seqs;
	vector<string> mlcs;
	vector< vector< vector<int> > > SucTabs;
	BiHashTable DM;
	vector<int> ID;
	vector<int> tlevel;
    vector< list<int> > precursor;
    vector<bool> inLCSPath;

};

int exe_rlpmlcs(vector<string>& seqs, string& alphasets, ostream& os, string& algo);

#endif // RLP_MLCS_H