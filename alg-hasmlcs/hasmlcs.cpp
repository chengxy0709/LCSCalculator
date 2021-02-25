#include "hasmlcs.h"

static int maxlevel = 0;
static bool opt = false;
static int cnt = 0;
static clock_t MAXTIME = 60 * 1000 * 1000;
static int optimal_len = INT32_MAX;

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

    int n = 0;
    for(int i = 0; i < seqs.size(); i++) if(seqs[i].length() > n) n = seqs[i].length();
    P = cal_pr_val(n, n, cmap.size());
	
}

void HASMLCS::run_for_BS(){

    priority_queue<Point<CordType>*, vector<Point<CordType>*>, priority1> Q;
    Point<CordType>* p0 = new Point<CordType>(seqs.size(), false, 0); 
    int loop = 0;

    maxlevel = 0;
    opt = false;
    PSETATTR(p0, new HASAttr(UB(p0), 0, 0, 0));
    ATTR(HASAttr, p0)->flag1 = 0;
    Q.push(p0);
    pset.insert(p0);
    while(Q.size() > 0 && !opt && MAXTIME > clock() /* not memory limit */){

        loop++;
        // init beam data
        vector< Point<CordType>* > B;

        for(int i = 0; i < beta && i < Q.size(); i++){
            Point<CordType> *p = pop_valid1(Q);
            if(p){
                ATTR(HASAttr, p)->flag1 = -1; // let it invalid in Q
                B.push_back(p);
            }
            else break;
        }

        // Beam Search
        while(!opt && B.size() > 0 && MAXTIME > clock()){
            HashSet Vext; 
            for(int i = 0; i < B.size(); i++){
                vector< Point<CordType>* > temp;
                ATTR( HASAttr, B[i] )->flag1 = -1;
                temp = ExpandNode(B[i], Q);
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
            Filter(Vext, Kfilter);
            /*for(auto p : Ref){
                p->print(seqs.size(), '\t');
                cout << "f:" << ATTR(HASAttr, p)->f << "g:" << ATTR(HASAttr, p)->g << "h:" << ATTR(HASAttr, p)->f - ATTR(HASAttr, p)->g << "\n";
            }
            cout << "===Ref\n";*/
            // Select beta points in Vext according to UB
            B = Reduce(Vext, beta);
            /*for(auto p : B){
                p->print(seqs.size(), '\t');
                cout << "f:" << ATTR(HASAttr, p)->f << "g:" << ATTR(HASAttr, p)->g << "h:" << ATTR(HASAttr, p)->f - ATTR(HASAttr, p)->g << "\n";
            }
            cout << "===B\n";*/
        }

        // A* Search
        int iter = 0;
        while(!opt && iter < delta && MAXTIME > clock() /*not memory limit */){
            Point<CordType> *point = pop_valid1(Q);
            if(point){
                ATTR(HASAttr, point)->flag1 = -1;
                ExpandNode(point, Q);
            }
            iter++;
        } 

    }
    
    cout << "The remain number of points in Q : " << Q.size() << endl;
    cout << "The number of points in set : " << pset.size() << endl;
    cout << "The number of expanding : " << cnt << endl;
    cout << "The number of loop : " << loop << endl;
}

void HASMLCS::run_for_ACS(){

    queue_type1 Q;
    vector< queue_type2 > Qlev;
    Point<CordType>* p0 = new Point<CordType>(seqs.size(), false, 0); 
    int loop = 0;
    int level = 0;

    maxlevel = 0;
    opt = false;
    PSETATTR(p0, new HASAttr(UB(p0), 0, 0, 0));
    ATTR( HASAttr, p0 )->flag1 = 0;
    ATTR( HASAttr, p0 )->flag2 = 0;
    Q.push(p0);
    Qlev.push_back(queue_type2{});
    Qlev[0].push(p0);
    pset.insert(p0);
    while(Q.size() > 0 && !opt && MAXTIME > clock() /* not memory limit */){

        loop++;
        level = 0;
        // perform ACS iteration
        while(!opt && level < Qlev.size() && MAXTIME > clock()){
            HashSet Vext;
            for(int i = 0; Qlev[level].size() > 0 && i < beta; i++){
                vector< Point<CordType>* > temp;
                Point<CordType> *point = pop_valid2(Qlev[level], level);
                if(point){
                    ATTR( HASAttr, point )->flag1 = -1;
                    ATTR( HASAttr, point )->flag2 = -1;
                    temp = ExpandNode(point, Q, Qlev);
                    // ExpandNode and Store children in Vext
                    Vext.insert(temp.begin(), temp.end());
                }
            }
            /*for(auto p : Vext){
                p->print(seqs.size(), '\t');
                cout << "f:" << ATTR(HASAttr, p)->f << "g:" << ATTR(HASAttr, p)->g << "h:" << ATTR(HASAttr, p)->f - ATTR(HASAttr, p)->g << "\n";
            }
            cout << "===Vext\n";*/
            // Select Kfilter most promising points as a reference set, and delete other points 
            // that can be dominated by these selected points.
            Filter(Vext, Kfilter);
            /*for(auto p : Vext){
                p->print(seqs.size(), '\t');
                cout << "f:" << ATTR(HASAttr, p)->f << "g:" << ATTR(HASAttr, p)->g << "h:" << ATTR(HASAttr, p)->f - ATTR(HASAttr, p)->g << "\n";
            }
            cout << "===B\n";*/
            level++;
        }

        // A* Search
        int iter = 0;
        while(!opt && iter < delta && MAXTIME > clock() /*not memory limit */){
            Point<CordType> *point = pop_valid1(Q);
            if(point){
                ATTR( HASAttr, point )->flag1 = -1;
                ATTR( HASAttr, point )->flag2 = -1;
                ExpandNode(point, Q);
            }
            iter++;
        } 

    }
    
    cout << "The remain number of points in Q : " << Q.size() << endl;
    cout << "The number of points in set : " << pset.size() << endl;
    cout << "The number of expanding : " << cnt << endl;
    cout << "The number of loop : " << loop << endl;
}

void HASMLCS::Filter(HashSet& Vext, int Kfilter){

    priority_queue< Point<CordType>*, vector< Point<CordType>* >, priority3 > Q;

    for(auto p : Vext){
        if(Q.size() > Kfilter){
            Q.pop();
        }
        Q.push(p);
    }

    while(Q.size() > Kfilter) Q.pop();
    while(Q.size() > 0){
        auto i = Vext.begin();
        while(i != Vext.end() ){
            if(is_successor(*i, Q.top(), seqs.size())){
                i = Vext.erase(i);
            }
            else i++;
        }
        Q.pop();
    }

}

bool UBComp(const Point<CordType> *p, const Point<CordType> *q){
    return (ATTR(HASAttr, p)->f - ATTR(HASAttr, p)->g) > (ATTR(HASAttr, q)->f - ATTR(HASAttr, q)->g);
}

vector< Point<CordType>* > HASMLCS::Reduce(HashSet& Vext, int beta){

    vector< Point<CordType>* > res;

    if(beta > Vext.size()){
        res = vector< Point<CordType>* >(Vext.begin(), Vext.end());
    }
    else{
        res = vector< Point<CordType>* >(Vext.begin(), Vext.end());
        sort(res.begin(), res.end(), UBComp);
        res.erase(res.begin() + beta, res.end());
    }

    return res;

}

vector< Point<CordType>* > HASMLCS::ExpandNode(Point<CordType>* p, queue_type1& Q){

    vector< Point<CordType>* > Vext;
    bool hasSuccessor = false;
    int l = ATTR(HASAttr, p)->g; // the level of p

    cnt++;
    for(int i = 0; i < SucTabs[0].size(); i++){
        Point<CordType>* suc = successor(p, SucTabs, i);
        if(suc){
            hasSuccessor = true;
            auto res = pset.insert(suc);
            if(res.second){ // new point
                int ub = UB(suc);
                PSETATTR(suc, (new HASAttr(ub + l + 1, l + 1, cal_k_norm(suc), 0, p)));
                ATTR(HASAttr, suc)->flag1 = 0;
                Q.push(suc); 
                Vext.push_back(suc);
            }
            else{ // update point
                Point<CordType>* ep = *res.first;
                if(ATTR(HASAttr, ep)->g < l + 1){ // a better solution            
                    ATTR(HASAttr, ep)->f += l + 1 - ATTR(HASAttr, ep)->g;
                    ATTR(HASAttr, ep)->g = l + 1;
                    ATTR(HASAttr, ep)->parent = p;
                    // Actually, here should update this point's position in Q.
                    Q.push(ep); 
                    ATTR(HASAttr, ep)->flag1 = 0; // let this point valid
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

    if(Q.size() <= 0 || maxlevel >= Qmax(Q)){
        opt = true;
    }
    if(maxlevel >= optimal_len) opt = true;
    
    return Vext;

}

vector< Point<CordType>* > HASMLCS::ExpandNode(Point<CordType>* p, queue_type1& Q, vector<queue_type2>& Qlev){

    vector< Point<CordType>* > Vext;
    bool hasSuccessor = false;
    int l = ATTR(HASAttr, p)->g; // the level of p

    cnt++;
    if(l + 1 >= Qlev.size()){
        Qlev.push_back(queue_type2{});
    }
    for(int i = 0; i < SucTabs[0].size(); i++){
        Point<CordType>* suc = successor(p, SucTabs, i);
        if(suc){
            hasSuccessor = true;
            auto res = pset.insert(suc);
            if(res.second){ // new point
                int ub = UB(suc);
                PSETATTR(suc, (new HASAttr(ub + l + 1, l + 1, cal_k_norm(suc), cal_ex_val(suc), p)));
                ATTR(HASAttr, suc)->flag1 = 0;
                ATTR(HASAttr, suc)->flag2 = l + 1;
                Q.push(suc); 
                Qlev[l + 1].push(suc);
                Vext.push_back(suc);
            }
            else{ // update point
                Point<CordType>* ep = *res.first;
                if(ATTR(HASAttr, ep)->g < l + 1){ // a better solution            
                    ATTR(HASAttr, ep)->f += l + 1 - ATTR(HASAttr, ep)->g;
                    ATTR(HASAttr, ep)->g = l + 1;
                    ATTR(HASAttr, ep)->parent = p;
                    // Actually, here should update this point's position in Q.
                    Q.push(ep); 
                    Qlev[l + 1].push(ep);
                    ATTR(HASAttr, ep)->flag1 = 0; // let this point valid
                    ATTR(HASAttr, ep)->flag2 = l + 1;
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

    if(Q.size() <= 0 || maxlevel >= Qmax(Q)){
        opt = true;
    }
    if(maxlevel >= optimal_len) opt = true;
    
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

vector< vector<double> > HASMLCS::cal_pr_val(int p, int q, int sigSize){

    vector< vector<double> > P = vector< vector<double> >(p, vector<double>(q, 0));
    for(int i = 0; i < p; i++){
        for(int j = 0; j < q; j++){
            if(i == 0) P[i][j] = 1;
            else if(i > j) P[i][j] = 0;
            else P[i][j] = (1.0 / sigSize) * P[i - 1][j - 1] + ((sigSize - 1) * 1.0 / sigSize) * P[i][j - 1];
        }
    }
    return P;

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

        cout << "The limit time (seconds) > ";
		cin >> MAXTIME;
        if(MAXTIME <= 0) MAXTIME = 60;
        MAXTIME *= 1000 * 1000;
		cout << "The best solution's length (set -1 for disabling this condition)> ";
		cin >> optimal_len;
        if(optimal_len < 0) optimal_len = INT32_MAX; 
		
        g_point_size = seqs.size();
		HASMLCS hasmlcs(seqs, alphasets, beta, delta, k);
		string lcs;
		clock_t start_t, end_t;
		
		start_t = clock();
        MAXTIME += start_t;
		hasmlcs.run_for_ACS();
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
