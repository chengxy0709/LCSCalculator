#include "quickdp.h"

QuickDP::QuickDP(vector<string> seqs, string alphabets)
 : seqs(seqs),alphaSize(alphabets.length())
{

	phash_init(seqs.size());
	cmap = build_alphabet_map(alphabets);
	SucTabs = cal_suc_tabs(seqs, cmap);
	
}

void QuickDP::run(){

    Point<CordType>* p0 = new Point<CordType>(seqs.size(), false, 0);
    int k = 0;
    vector< HashSet > D;

    D.push_back(HashSet{p0});

    while(D[k].size() != 0){
        vector< HashSet > Pars = vector< HashSet >(alphaSize, HashSet());
        HashSet trash;
        D.push_back(HashSet());
        for(Point<CordType> *point : D[k]){
            vector<Point<CordType>*> Par;
            // generate successors
            for(int i = 0; i < alphaSize; i++){
                Point<CordType> *suc = successor(point, SucTabs, i);
                if(suc){
                    if(trash.find(suc) == trash.end())
                        Par.push_back(suc);
                    else delete suc;
                }
            }
            // minimal Par
            vector<Point<CordType>*> minPar = Minima(Par, seqs.size());
            // put all the points belong to minPar in Pars correspondingly, meanwhile, put all the points
            // belong to (Par - minPar) in trash.
            sort(minPar.begin(), minPar.end());
            sort(Par.begin(), Par.end());
            for(int i = 0, j = 0; i < Par.size(); i++){
                if(j < minPar.size() && Par[i] == minPar[j]){
                    auto res = Pars[CharIndex(minPar[j])].insert(minPar[j]);
                    if(!res.second) delete minPar[j];
                    j++;
                }
                else{
                    int idx = CharIndex(Par[i]);
                    auto iter = Pars[idx].find(Par[i]);
                    if(iter != Pars[idx].end()) {
                        delete *iter;
                        Pars[idx].erase(iter);
                    }
                    auto res = trash.insert(Par[i]);
                    if(!res.second) delete Par[i];
                }
            }
        }
        // calculate Minimal Pars
        for(int i = 0; i < alphaSize; i++){
            vector<Point<CordType>*> Pars_vec = set2vec(Pars[i]);
            Pars_vec = Minima(Pars_vec, seqs.size());
            for(auto p : Pars_vec) D[k + 1].insert(p);
            for(auto p : Pars[i]){ // delete dominated point
                if(D[k + 1].find(p) == D[k + 1].end()) delete p;
            }
        }
        k++;
    }

    //get a lcs
    k--;
    Point<CordType> *curPoint = *D[k].begin();
    char c = seqs[0][curPoint->cord[0] - 1];
    lcs = "";
	lcs += c;
    while(k > 1){
        for(auto point : D[k - 1]){
            Point<CordType> *suc = successor(point, SucTabs, cmap[c]);
            if(suc == nullptr) continue;
            if(Key_equal()(suc, curPoint)){
                curPoint = point;
                c = seqs[0][curPoint->cord[0] - 1];
                delete suc; break;
            }
            else delete suc;
        }
        lcs = c + lcs;
        k--;
    }
	
	// free memory
	for(auto s : D){
		for(auto p : s) delete p;
	}

}

vector< Point<CordType>* > QuickDP::Minima(vector< Point<CordType>* >& points, int dim){

    Qsort(points, 0, points.size() - 1, 0); // sorted by x-axis
    return Divide(points, dim);

}

inline bool QuickDP::analyzePoints(vector< Point<CordType>* >& points, int dim){
    int refer = points[0]->cord[dim];
    for(auto point : points){
        if(point->cord[dim] != refer) return true;
    }
    return false; // has no relation of domination.
}

vector<Point<CordType>*> QuickDP::Divide(vector< Point<CordType>* >& points, int d){

    vector<Point<CordType>*> res;

    // only a or zero point
    if ( points.size() <= 1 ){
        return points;
    }

    if (d == 2){ // conquer
        int ymin = INT_MAX;
        for(auto point : points){
            if(point->cord[1] < ymin){
                ymin = point->cord[1];
                if(res.size() > 0 && point->cord[0] == res.back()->cord[0])
                    res.pop_back();
                res.push_back(point);
            }
        }

        return res;
    }

    if(!analyzePoints(points, d - 1)){
        return Divide(points, d - 1);
    }

    // divide
    vector<Point<CordType>*> A;
    vector<Point<CordType>*> B;
    // find median number
    int median = Qmedian(points, d - 1);
    int maxval = vmax(points, d - 1);
    // record the index of every point and place it in A or B
    bool flag = (median == maxval) ? true : false;
    for(auto point : points){
        if(point->cord[d - 1] < median){
            A.push_back(point);
        }
        else if(point->cord[d - 1] > median){
            B.push_back(point);
        }
        else{
            if(flag) B.push_back(point);
            else A.push_back(point);
        }
    }

    vector<Point<CordType>*> Ares = Divide(A, d);
    vector<Point<CordType>*> Bres = Divide(B, d);

    // merge
    // set label
    for(auto p : Ares){
		SETATTRINT(p, 1);
    }
    for(auto p : Bres){
		SETATTRINT(p, 0);
    }
    Ares = mergeSortedVecter(Ares, Bres, 0);
    // delete the point in Ares dominated by any point in Bres
    res = Union(Ares, d - 1);

    return res;

}

vector<Point<CordType>*> QuickDP::Union(vector< Point<CordType>* >& points, int d){

    vector<Point<CordType>*> res;

    // only a or zero point
    if ( points.size() <= 1){
        return points;
    }

    if (d == 2){ // conquer
        int ymin = INT_MAX;
        for(auto point : points){
            if(ATTRINT(point) != 0){ // label is A
                if(point->cord[1] < ymin){
                    ymin = point->cord[1];
                    if(res.size() > 0 &&
                       point->cord[0] == res.back()->cord[0] &&
                       point->cord[1] <= res.back()->cord[1] &&
                       ATTRINT(res.back()) == 0){
                        res.pop_back();
                    }
                }
                res.push_back(point);
            }
            else{
                if(point->cord[1] < ymin){
                    res.push_back(point);
                }
            }
        }
        return res;
    }

    if(!analyzePoints(points, d - 1)){
        return Union(points, d - 1);
    }

    // divide
    vector<Point<CordType>*> A;
    vector<Point<CordType>*> B;
    // find median number
    int median = Qmedian(points, d - 1);
    int maxval = vmax(points, d - 1);
    // record the index of every point and place it in A or B
    bool flag = (median == maxval) ? true : false;
    for(auto point : points){
        if(point->cord[d - 1] < median){
            A.push_back(point);
        }
        else if(point->cord[d - 1] > median){
            B.push_back(point);
        }
        else{
            if(flag) B.push_back(point);
            else A.push_back(point);
        }
    }

    vector<Point<CordType>*> Ares = Union(A, d);
    vector<Point<CordType>*> Bres = Union(B, d);

    // delete the point labeled 0 in Bres dominated by any point labeled 1 in Ares
    vector<Point<CordType>*> UA;
    vector<Point<CordType>*> UB;
    vector<Point<CordType>*> Blabeled1; // the point labeled 1 in Bres
    for(auto p : Bres){
        if(ATTRINT(p) == 0) UB.push_back(p);
        else Blabeled1.push_back(p);
    }
    for(auto p : Ares){
        if(ATTRINT(p) == 1) UA.push_back(p);
    }
    UA = mergeSortedVecter(UA, UB, 0);
    vector<Point<CordType>*> Ures = Union(UA, d - 1);

    // merge consequence
    Ares = mergeSortedSet(Ares, Blabeled1, 0, seqs.size());
    Ares = mergeSortedSet(Ares, Ures, 0, seqs.size());

    return Ares;

}

int exe_quickdp(vector<string>& seqs, string& alphasets, ostream& os, string& algo){
    if(algo == QUICKDPSYM){
        QuickDP quickdp(seqs, alphasets);
		string lcs;
		clock_t start_t, end_t;
		
		start_t = clock();
		quickdp.run();
		end_t = clock();
		lcs = quickdp.LCS();
		os << "Result(by " << algo << "):\n";
		os << "time(us) : " << end_t - start_t << "\n";
		os << "the length of lcs : " << lcs.length() << "\n";
		os << "a lcs : " << lcs << "\n";
        return 0;
    }
    else{
        return -1;
    }
}