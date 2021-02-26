#ifndef TOP_MLCS_H
#define TOP_MLCS_H

#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <list>
#include <algorithm>
#include <time.h>

#include "point.h"
#include "phash.h"
#include "tool.h"

#define ALG_REQUIRED

#define TOPMLCSSYM "TOPMLCS"

class TOP_MLCS{

public:    
    TOP_MLCS(vector<string>& seqs, string& alphabets);
    ~TOP_MLCS();
    
    void run();
    vector<string> MLCS() const {return mlcs;}

private:
    int construct_ICSG();
    int ForwardTopSort(int maxIndex);
    void BackwardTopSort(int maxlevel, int pinf_index);
    void GetMLCS(string& LCSRecord, int index);

    vector<string> seqs;
    vector<string> mlcs;
    vector< vector< vector<int> > > SucTabs;
    BiHashTable DM;
    vector<int> ID;
    vector<int> tlevel;
    vector< list<int> > precursor;

};

int exe_topmlcs(vector<string>& seqs, string& alphasets, ostream& os, string& algo);
void UsageforTOPMLCS();

#endif // TOP_MLCS_H