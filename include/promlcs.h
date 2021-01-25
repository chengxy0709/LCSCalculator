#ifndef PRO_MLCS_H
#define PRO_MLCS_H

#include <iostream>
#include <list>
#include <string>
#include <unordered_map>
#include <queue>
#include <limits.h>

#include "point.h"
#include "tool.h"
#include "dtree.h"

#define ALG_REQUIRED

#define PROMLCSSYM "PROMLCS"

using namespace std;

class PRO_MLCS{
	
public:
	PRO_MLCS(vector<string> seqs, string alphabets, int scnum = -1);
	~PRO_MLCS(){}
	
	void run();
	string LCS() const{return lcs;}
	
private:
	void calSU(unordered_map<char, int>& cmap);
	int calg(Point<CordType> *p, int order);
	int caldist(Point<CordType> *p);

	vector<string> seqs;
	int alphaSize;
	vector< vector< vector<int> > > SucTabs;
	vector< vector< vector<int> > > SU; // for estimating the upper bound of a local LCS
	string lcs;
	
	int scnum; // the number of points that calculate in a single iteration
	
};

int exe_promlcs(vector<string>& seqs, string& alphasets, ostream& os, string& algo);

#endif //PRO_MLCS_H