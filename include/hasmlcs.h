#ifndef HYBRID_A_STAR_H
#define HYBRID_A_STAR_H

#include <iostream>
#include <algorithm>
#include <queue>
#include <list>
#include <time.h>

#include "point.h"
#include "phash.h"
#include "tool.h"

#define ALG_REQUIRED

#define HASMLCSSYM "HASMLCS"

using namespace std;


// priority definition
struct priority1{
    bool comp(const Point<CordType>* p1, const Point<CordType>* p2) const;
    bool operator () (const Point<CordType>* p1, const Point<CordType>* p2) const{
        return comp(p1, p2);
    }
};
struct priority2{
    bool comp(const Point<CordType>* p1, const Point<CordType>* p2) const;
    bool operator () (const Point<CordType>* p1, const Point<CordType>* p2) const{
        return comp(p1, p2);
    }
};
struct priority3{
    bool comp(const Point<CordType>* p1, const Point<CordType>* p2) const;
    bool operator () (const Point<CordType>* p1, const Point<CordType>* p2) const{
        return comp(p1, p2);
    }
};

struct HASAttr{
	int f;
	int g;
    double k; // k-norm
    double EX; // expected length
	Point<CordType>* parent;
	int flag1;
	int flag2;
	HASAttr() : f(0),g(0),k(0),EX(0),parent(NULL){}
	HASAttr(int f, int g, double k, double EX, Point<CordType>* parent = NULL) 
        : f(f),g(g),k(k),EX(EX),parent(parent){}
};

class HASMLCS{

public:
	HASMLCS(vector<string>& seqs, string& alphabets, int beta = -1, int delta = -1, int K = -1);
	~HASMLCS(){

	}
	
	void run_for_BS();
	void run_for_ACS();
	string LCS() const {return lcs;}

private:

	vector< vector<double> > cal_pr_val(int p, int q, int sigSize);
    double cal_k_norm(const Point<CordType>* p);
    double cal_ex_val(const Point<CordType>* p);
	int UB(Point<CordType>* p);
	void Filter(HashSet& Vext, int Kfilter);
	vector< Point<CordType>* > Reduce(HashSet& Vext, int beta);
	int getLCS( Point<CordType>* p );

	template<class T>
	Point<CordType>* pop_valid1(T& Q){
		Point<CordType> *p = NULL;
		while(Q.size() > 0) {
			p = Q.top();
			Q.pop();
			if(ATTR(HASAttr, p)->flag1 == 0) break;
		} 
		return (ATTR(HASAttr, p)->flag1 == 0) ? p : NULL;
	}
	template<class T>
	Point<CordType>* pop_valid2(T& Q, int level){
		Point<CordType> *p = NULL;
		while(Q.size() > 0) {
			p = Q.top();
			Q.pop();
			if(ATTR(HASAttr, p)->flag2 == level) break;
		} 
		return (ATTR(HASAttr, p)->flag2 == level) ? p : NULL;
	}
	template<class T>
	int Qmax(T& Q){
        Point<CordType> *p = NULL;
		while(Q.size() > 0) {
			p = Q.top();
			if(ATTR(HASAttr, p)->flag1 == 0) break;
			Q.pop();
		}
		return (ATTR(HASAttr, p)->flag1 == 0) ? ATTR(HASAttr, p)->f : 0;
	}

	using queue_type1 = priority_queue<Point<CordType>*, vector<Point<CordType>*>, priority1>;
	using queue_type2 = priority_queue<Point<CordType>*, vector<Point<CordType>*>, priority2>;
	using queue_type3 = priority_queue<Point<CordType>*, vector<Point<CordType>*>, priority3>;

	vector< Point<CordType>* > ExpandNode(Point<CordType>* p, queue_type1& Q);
	vector< Point<CordType>* > ExpandNode(Point<CordType>* p, queue_type1& Q, vector<queue_type2>& Qlev);

	vector<string> seqs;
	vector< vector< vector<int> > > SucTabs;
	vector< vector< vector<int> > > ScoreTabs; // for estimating the length of a local LCS
	vector< vector< vector<int> > > CountTabs; // for estimating the length of a local LCS
	HashSet pset;
	vector< vector<double> > P; // pr table

	string lcs;

	int beta; // BS width
	int delta; // the number of iteration of A* Search
	int Kfilter;
	
};

int exe_hasmlcs(vector<string>& seqs, string& alphasets, ostream& os, string& algo);

#endif // HYBRID_A_STAR_H
