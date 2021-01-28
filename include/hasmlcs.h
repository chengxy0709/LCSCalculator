#ifndef HYBRID_A_STAR_H
#define HYBRID_A_STAR_H

#include <iostream>
#include <algorithm>
#include <queue>
#include <list>
#include <time.h>
#include <set>

#include "point.h"
#include "phash.h"
#include "tool.h"

#define ALG_REQUIRED

#define HASMLCSSYM "HASMLCS"

using namespace std;

struct HASAttr{
	int f;
	int g;
    double k; // k-norm
    double EX; // expected length
	Point<CordType>* parent;
	HASAttr() : f(0),g(0),k(0),EX(0),parent(NULL){}
	HASAttr(int f, int g, double k, double EX, Point<CordType>* parent = NULL) 
        : f(f),g(g),k(k),EX(EX),parent(parent){}
};


// priority definition
struct priority1{
    bool comp(const Point<CordType>* p1, const Point<CordType>* p2) const;
    bool operator () (const Point<CordType>* p1, const Point<CordType>* p2) const{
        return comp(p1, p2);
    }
};
struct priority2{
    bool comp(const Point<CordType>* p1, const Point<CordType>* p2) const;
    bool operator () (const Point<CordType>* p1, const Point<CordType>* p2) const;
};
struct priority3{
    bool comp(const Point<CordType>* p1, const Point<CordType>* p2) const;
    bool operator () (const Point<CordType>* p1, const Point<CordType>* p2) const{
        return comp(p1, p2);
    }
};

class HASMLCS{

public:
	HASMLCS(vector<string>& seqs, string& alphabets, int beta = -1, int delta = -1, int K = -1);
	~HASMLCS(){

	}
	
	void run();
	string LCS() const {return lcs;}

private:

    double cal_k_norm(const Point<CordType>* p);
    double cal_ex_val(const Point<CordType>* p);
	int UB(Point<CordType>* p);
	vector< Point<CordType>* > Filter(HashSet& Vext, int Kfilter);
	vector< Point<CordType>* > Reduce(HashSet& Vext, vector< Point<CordType>* >& Ref, int beta);
	int getLCS( Point<CordType>* p );

	vector< Point<CordType>* > ExpandNode(Point<CordType>* p, set<Point<CordType>*, priority3>& Q);

	vector<string> seqs;
	vector< vector< vector<int> > > SucTabs;
	vector< vector< vector<int> > > ScoreTabs; // for estimating the length of a local LCS
	vector< vector< vector<int> > > CountTabs; // for estimating the length of a local LCS
	HashSet pset;

	string lcs;

	int beta; // BS width
	int delta; // the number of iteration of A* Search
	int Kfilter;
	
};

int exe_hasmlcs(vector<string>& seqs, string& alphasets, ostream& os, string& algo);

#endif // HYBRID_A_STAR_H