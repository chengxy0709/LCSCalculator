#ifndef TOOL_H
#define TOOL_H

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>

#include "point.h"

using namespace std;

unordered_map<char, int> build_alphabet_map(string& alphabets);

vector< vector<int> > cal_suc_tab(string& seq, unordered_map<char, int>& cmap, bool disp = false);
vector< vector< vector<int> > > cal_suc_tabs(vector<string>& seqs, unordered_map<char, int>& cmap, bool disp = false);

vector< vector<int> > cal_count_tab(string& seq, unordered_map<char, int>& cmap, bool disp = false);
vector< vector< vector<int> > > cal_count_tabs(vector<string>& seqs, unordered_map<char, int>& cmap, bool disp = false);
int UpperBound_by_CountTabs(Point<CordType> *p, vector< vector< vector<int> > >& CountTabs);

vector< vector< vector<int> > > cal_score_tabs(vector<string>& seqs);
int UpperBound_by_ScoreTabs(Point<CordType> *p, vector< vector< vector<int> > >& ScoreTabs);

Point<CordType>* successor(Point<CordType>* p, vector< vector< vector<int> > >& SucTabs, int i);
bool is_successor(Point<CordType>* p, Point<CordType>* q, int num);
bool is_immediate_successor(Point<CordType>* p, Point<CordType>* q , vector< vector< vector<int> > >& SucTabs);

void Qsort(vector< Point<CordType>* >& arr, int low, int high, int dim);
int Qmedian(vector< Point<CordType>* >& arr, int dim);
int vmax(vector< Point<CordType>* >& arr, int dim);

// A and B are two sorted vectors in dim, and the function mergeSortedVecter can merge the two vectors
// to A by this order.
vector<Point<CordType>*> mergeSortedVecter(vector<Point<CordType>*>& A, vector<Point<CordType>*>& B, int dim);
vector<Point<CordType>*> mergeSortedSet(vector<Point<CordType>*>& A, vector<Point<CordType>*>& B, int dim, int pointsize);

#endif //TOOL_H