#ifndef QUICK_DP_H
#define QUICK_DP_H

#include <iostream>
#include <vector>
#include <limits.h>

#include "point.h"
#include "tool.h"
#include "phash.h"

#define ALG_REQUIRED

#define QUICKDPSYM "QUICKDP"

class QuickDP{
	
public:
	QuickDP(vector<string> seqs, string alphabets);
	~QuickDP(){}
	
	vector< Point<CordType>* > Minima(vector< Point<CordType>* >& points, int dim);	
	void run();
	string LCS() const {return lcs;}
	
private:
	inline bool analyzePoints(vector< Point<CordType>* >& points, int dim);
	vector<Point<CordType>*> Divide(vector< Point<CordType>* >& points, int d);
	vector<Point<CordType>*> Union(vector< Point<CordType>* >& points, int d);
	inline int CharIndex(const Point<CordType> *p){return cmap[seqs[0][p->cord[0] - 1]];}

	vector<string> seqs;
	unordered_map<char, int> cmap;
	int alphaSize;
	vector< vector< vector<int> > > SucTabs;
	string lcs;
	
};

int exe_quickdp(vector<string>& seqs, string& alphasets, ostream& os, string& algo);

#endif //QUICK_DP_H