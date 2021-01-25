#ifndef MLCSAPP_H
#define MLCSAPP_H

#include <iostream>
#include <algorithm>

#include "point.h"
#include "phash.h"
#include "tool.h"

#define ALG_REQUIRED

#define MLCSAPPSYM "MLCSAPP"

using namespace std;

struct Attribute{
	int f;
	int g;
	Point<CordType>* parent;
	Attribute() : f(0),g(0),parent(NULL){}
	Attribute(int f, int g, Point<CordType>* parent = NULL) : f(f),g(g),parent(parent){}
};

class MLCSAPP{

public:
	MLCSAPP(vector<string>& seqs, string& alphabets, int k = -1, int c = -1);
	~MLCSAPP();
	
	void run();
	string LCS() const {return lcs;}
private:
	void compute_score_tabs();
	void computeF(Point<CordType> *p);
	int getLCS( Point<CordType>* p );

	vector<string> seqs;
	vector< vector< vector<int> > > SucTabs;
	vector< vector< vector<int> > > ScoreTabs; // for estimating the length of a local LCS
	string lcs;
	
	int CONSTK;
	int CONSTC;
};

int exe_mlcsapp(vector<string>& seqs, string& alphasets, ostream& os, string& algo);

#endif // MLCSAPP_H