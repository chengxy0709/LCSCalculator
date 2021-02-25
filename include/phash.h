#ifndef PHASH_H
#define PHASH_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdlib.h>
#include <math.h>

#include "point.h"

using namespace std;
void phash_init(int seqnum, int fac = 10000000);
struct Key_equal{
    bool operator() (const Point<CordType> *p1, const Point<CordType> *p2) const;
};
struct Hash_Function{
    size_t operator() (const Point<CordType> *p) const;
};

using HashMap = unordered_map<Point<CordType>*, int, Hash_Function, Key_equal>;
using HashSet = unordered_set<Point<CordType>*, Hash_Function, Key_equal>;

vector<Point<CordType>*> set2vec(HashSet& pset);

class BiHashTable{

public:
    BiHashTable(int volume = -1);
    ~BiHashTable();

    bool insert(Point<CordType> *p);
    int at(Point<CordType> *p) {return Point2IndexTab[p];}
    Point<CordType>* at(int i) {return Index2PointTab[i];}
    unsigned int size() const{return sz;}
    unsigned int bucket_count() const{return Point2IndexTab.bucket_count();}
    unsigned int max_bucket_count() const{return Point2IndexTab.max_bucket_count();}
    unsigned int maxnum_in_each_bucket() const;

private:
    HashMap Point2IndexTab;
    vector< Point<CordType>* > Index2PointTab;
       unsigned int sz;     

};

#endif //PHASH_H
