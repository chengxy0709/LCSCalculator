#include "hasmlcs.h"

static int maxlevel = 0;
static int maxval = 0; // the max value f of a point in Q
static bool opt = false;

HASMLCS::HASMLCS(vector<string>& seqs, string& alphabets, int beta, int delta, int K){
	
	unordered_map<char, int> cmap = build_alphabet_map(alphabets);
	phash_init(seqs.size(), 1000000);
	this->seqs = seqs;
	if(beta >= 0) this->beta = beta;
	else this->beta = 1;
	if(delta >= 0) this->delta = delta;
	else this->delta = 1;
    if(K >= 0) Kfilter = K;
	else Kfilter = 0;
	SucTabs = cal_suc_tabs(seqs, cmap);
	ScoreTabs = cal_score_tabs(seqs);
    CountTabs = cal_count_tabs(seqs, cmap, false);
	
}

void HASMLCS::run(){

    set< Point<CordType>*,  priority3 > Q;
    Point<CordType>* p0 = new Point<CordType>(seqs.size(), false, 0); 

    maxlevel = 0;
    maxval = 0;
    opt = false;
    PSETATTR(p0, new HASAttr(UB(p0), 0, 0, 0));
    Q.insert(p0);
    pset.insert(p0);
    while(Q.size() > 0 && !opt /* not time limit and not memory limit */){

        // init beam data
        vector< Point<CordType>* > B;
        auto endp = Q.begin();

        maxval = ATTR(HASAttr, (*Q.begin()))->f;
        for(int i = 0; i < beta && i < Q.size(); i++) endp++;
        B.insert(B.end(), Q.begin(), endp);
        Q.erase(Q.begin(), endp);
        

        // Beam Search
        while(B.size() > 0){
            HashSet Vext; 
            for(int i = 0; i < B.size(); i++){
                vector< Point<CordType>* > temp = ExpandNode(B[i], Q);
                // ExpandNode and Store children in Vext
                Vext.insert(temp.begin(), temp.end());
            }
            /*for(auto p : Vext){
                p->print(seqs.size(), '\t');
                cout << "f:" << ATTR(HASAttr, p)->f << "g:" << ATTR(HASAttr, p)->g << "h:" << ATTR(HASAttr, p)->f - ATTR(HASAttr, p)->g << "\n";
            }
            cout << "===Vext\n";*/
            // Select Kfilter most promising points as a reference set, and delete other points 
            // that can be dominated by these selected points.
            vector< Point<CordType>* > Ref = Filter(Vext, Kfilter);
            /*for(auto p : Ref){
                p->print(seqs.size(), '\t');
                cout << "f:" << ATTR(HASAttr, p)->f << "g:" << ATTR(HASAttr, p)->g << "h:" << ATTR(HASAttr, p)->f - ATTR(HASAttr, p)->g << "\n";
            }
            cout << "===Ref\n";*/
            // Select beta points in Vext according to UB
            B = Reduce(Vext, Ref, beta);
            /*for(auto p : B){
                p->print(seqs.size(), '\t');
                cout << "f:" << ATTR(HASAttr, p)->f << "g:" << ATTR(HASAttr, p)->g << "h:" << ATTR(HASAttr, p)->f - ATTR(HASAttr, p)->g << "\n";
            }
            cout << "===B\n";*/
        }

        // A* Search
        int iter = 0;
        while(Q.size() > 0 && iter < delta /* not time limit and not memory limit */){
            Point<CordType> *point = *(Q.begin());
            Q.erase(Q.begin());
            ExpandNode(point, Q);
            iter++;
        }

    }

}

vector< Point<CordType>* > HASMLCS::Filter(HashSet& Vext, int Kfilter){

    priority_queue< Point<CordType>*, vector< Point<CordType>* >, priority3 > Q;
    vector< Point<CordType>* > res;

    for(auto p : Vext){
        if(Q.size() > Kfilter){
            Q.pop();
        }
        Q.push(p);
    }

    while(Q.size() > Kfilter) Q.pop();
    while(Q.size() > 0){
        res.push_back(Q.top());
        Q.pop();
    }

    return res;

}

bool UBComp(const Point<CordType> *p, const Point<CordType> *q){
    return (ATTR(HASAttr, p)->f - ATTR(HASAttr, p)->g) > (ATTR(HASAttr, q)->f - ATTR(HASAttr, q)->g);
}

vector< Point<CordType>* > HASMLCS::Reduce(HashSet& Vext, vector< Point<CordType>* >& Ref, int beta){

    vector< Point<CordType>* > res;
    vector< Point<CordType>* > temp(Vext.begin(), Vext.end());

    sort(temp.begin(), temp.end(), UBComp);
    for(int i = 0; i < temp.size(); i++){
        int j = 0;
        for(j = 0; j < Ref.size() ; j++){
            if(is_successor(temp[i], Ref[j], seqs.size())) break;
        }
        if(j >= Ref.size()) res.push_back(temp[i]);
        if(res.size() >= beta) break;
    }

    return res;

}

vector< Point<CordType>* > HASMLCS::ExpandNode(Point<CordType>* p, set<Point<CordType>*,  priority3>& Q){

    vector< Point<CordType>* > Vext;
    bool hasSuccessor = false;
    int l = ATTR(HASAttr, p)->g; // the level of p

    for(int i = 0; i < SucTabs[0].size(); i++){
        Point<CordType>* suc = successor(p, SucTabs, i);
        if(suc){
            hasSuccessor = true;
            auto res = pset.insert(suc);
            if(res.second){ // new point
                int ub = UB(suc);
                PSETATTR(suc, (new HASAttr(ub + l + 1, l + 1, cal_k_norm(suc), 0, p)));
                Q.insert(suc);
                Vext.push_back(suc);
            }
            else{ // update point
                Point<CordType>* ep = *(res.first);
                if(ATTR(HASAttr, ep)->g < l + 1){ // a better solution            
                    auto num = Q.erase(ep);
                    ATTR(HASAttr, ep)->f += l + 1 - ATTR(HASAttr, ep)->g;
                    ATTR(HASAttr, ep)->g = l + 1;
                    ATTR(HASAttr, ep)->parent = p;
                    // Actually, here should update this point's position in Q.
                    if(num > 0) Q.insert(ep);
                    if(num > 1) cout << num << endl;
                }
                Vext.push_back(ep);

                delete ATTR(HASAttr, suc);
                delete suc;
            }
        }
    }

    if(!hasSuccessor){ // p is a completed point
        // update the best solution
        if(ATTR(HASAttr, p)->g > maxlevel){
            maxlevel = getLCS(p);
        }
    }

    if(Q.size() <= 0 || maxlevel > ATTR(HASAttr, (*(Q.begin())))->f){
        opt = true;
    }
    
    return Vext;

}

int HASMLCS::UB(Point<CordType>* p){

    int UB1 = UpperBound_by_ScoreTabs(p, ScoreTabs);
    int UB2 = UpperBound_by_CountTabs(p, CountTabs);

    return min(UB1, UB2);

}

int HASMLCS::getLCS( Point<CordType>* p ){
    lcs.clear();
	while(p->cord[0] != 0){
		lcs = seqs[0][p->cord[0] - 1] + lcs;
		p = ATTR(HASAttr, p)->parent;
	}
	return lcs.length();
}

int exe_hasmlcs(vector<string>& seqs, string& alphasets, ostream& os, string& algo){
    if(algo == HASMLCSSYM){
		/* else argument*/
		int beta, delta, k;
		cout << "The const value beta(negative value respresents default value) > ";
		cin >> beta;
		cout << "The const value delta(negative value respresents default value) > ";
		cin >> delta;
		cout << "The const value Kfilter(negative value respresents default value) > ";
		cin >> k;
		
		HASMLCS hasmlcs(seqs, alphasets, beta, delta, k);
		string lcs;
		clock_t start_t, end_t;
		
		start_t = clock();
		hasmlcs.run();
		end_t = clock();
		lcs = hasmlcs.LCS();
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