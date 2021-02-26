#include <sstream>

#include "mlcsapp.h"

MLCSAPP::MLCSAPP(vector<string>& seqs, string& alphabets, int k, int c){
    
    unordered_map<char, int> cmap = build_alphabet_map(alphabets);
    phash_init(seqs.size(), 1000000);
    this->seqs = seqs;
    if(c > 0) CONSTC = c;
    else CONSTC = 20;
    if(k > 0) CONSTK = k;
    else CONSTK = 2000;
    SucTabs = cal_suc_tabs(seqs, cmap);
    compute_score_tabs();
    
}

void MLCSAPP::compute_score_tabs(){
    
    for (int i = 0; i < (seqs.size() + 1) / 2; i++) {
        int i1 = i << 1;
        int i2 = 1+ (i << 1);
        if (i2 >= seqs.size()) i2 = 0;

        ScoreTabs.push_back(vector< vector<int> >(seqs[i1].length() + 1, vector<int>(seqs[i2].length() + 1, 0)));
        for (int x = seqs[i1].length() - 1; x >= 0; x--) {
            for (int y = seqs[i2].length() - 1; y >= 0; y--) {
                if (seqs[i1][x] == seqs[i2][y])
                    ScoreTabs[i][x][y] = ScoreTabs[i][x + 1][y + 1] + 1;
                else
                    ScoreTabs[i][x][y] = max(ScoreTabs[i][x+1][y], ScoreTabs[i][x][y+1]);
            }
        }

    }
    
}

void MLCSAPP::computeF(Point<CordType> *p){

    int h = ScoreTabs[0][0][0];
    for (int i = 0; i < (seqs.size() + 1) / 2; i++) {
        int i1 = i << 1;
        int i2 = 1+ (i << 1);
        if (i2 >= seqs.size()) i2 = 0;
        h = min(h, ScoreTabs[i][p->cord[i1]][p->cord[i2]]);
    }

    ATTR(Attribute, p)->f = h + ATTR(Attribute, p)->g;
    
}


int MLCSAPP::getLCS( Point<CordType>* p ) {
    while(p->cord[0] != 0){
        lcs = seqs[0][p->cord[0] - 1] + lcs;
        p = ATTR(Attribute, p)->parent;
    }
    return lcs.length();
}

struct Greater_Fun{
    size_t operator() (const Point<CordType>* a, const Point<CordType>* b) const{
        return ATTR(Attribute, a)->f > ATTR(Attribute, b)->f;
    }
};

void MLCSAPP::run(){

    HashSet nodeset(CONSTK);
    vector< vector< Point<CordType>* > > trash;
    Point<CordType> *p0 = new Point<CordType>(seqs.size(), false, 0);
    bool flag = true;
    
    PSETATTR(p0, new Attribute);
    computeF(p0);
    nodeset.insert(p0);

    while(nodeset.size() > 0) {
        vector< Point<CordType>* > nodeset_t;
        int maxf = 0;

        for(auto& p : nodeset){
            if(ATTR(Attribute, p)->f > maxf) maxf = ATTR(Attribute, p)->f;
        }
        for(auto& p : nodeset){
            if(ATTR(Attribute, p)->f >= (maxf - CONSTC)){
                nodeset_t.push_back(p);
            }
            else{
                delete ATTR(Attribute, p);
                delete p;
            }
        }

        if(nodeset_t.size() > CONSTK){
            sort(nodeset_t.begin(), nodeset_t.end(), Greater_Fun());
            for(unsigned int i = nodeset_t.size() - 1; i > CONSTK; i--){
                delete ATTR(Attribute, nodeset_t.back());
                delete nodeset_t.back();
                nodeset_t.pop_back();
            }
        }

        trash.push_back(nodeset_t);
        nodeset.clear();
        flag = true;
        for(auto& p : nodeset_t){
            if (ATTR(Attribute, p)->f - ATTR(Attribute, p)->g != 0) { // h(p) != 0
                flag = false;
            }

            for(int i = 0; i < SucTabs[0].size(); i++) {
                Point<CordType> *suc = successor(p, SucTabs, i);
                if(suc){
                    auto q = nodeset.find(suc);
                    if (q != nodeset.end()) { //node is already in node set
                        Attribute *sucAttr = ATTR(Attribute, (*q));
                        if ( sucAttr->g < ATTR(Attribute, p)->g + 1) {
                            sucAttr->parent = p;
                            sucAttr->g = ATTR(Attribute, p)->g + 1;
                            computeF(*q);
                        }
                        delete suc;
                    } 
                    else {
                        PSETATTR(suc, new Attribute(0, ATTR(Attribute, p)->g + 1, p));
                        computeF(suc);
                        nodeset.insert(suc);
                    }
                }
            }
        }
        if(flag || nodeset.size() == 0) break;
    }
    
    getLCS(trash.back()[0]);
    
    for(auto& vec : trash){
        for(auto& p : vec){
            delete ATTR(Attribute, p);
        }
    }
    
}

MLCSAPP::~MLCSAPP(){
        
}

int getparams(string& params, int& k, int& c){

    istringstream is(params);
    string opt;

    while(is >> opt){
        if(opt == "k") is >> k;
        else if(opt == "c") is >> c;
        else return -1;
    }

    return 0;

}

int exe_mlcsapp(vector<string>& seqs, string& alphasets, ostream& os, string& algo, string params){
    if(algo == MLCSAPPSYM){
        /* else argument*/
        int CONSTK, CONSTC;
        if(getparams(params, CONSTK, CONSTC)){
            cout << "extra parameter error" << endl;
            return 0;
        }
        
        MLCSAPP mlcsapp(seqs, alphasets, CONSTK, CONSTC);
        string lcs;
        clock_t start_t, end_t;
        
        start_t = clock();
        mlcsapp.run();
        end_t = clock();
        lcs = mlcsapp.LCS();
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

void UsageforMLCSAPP(){
    cout << endl;
    cout << "Information for MLCSAPP:\n" << endl;
    cout << "Description:" << endl;
    cout << "\tThis algorithm is re-implemented according to the article  \"A Fast Heuristic Search Algorithm for Finding \
the Longest Common Subsequence of Multiple Strings\". This is a variant for A* algorithm." << endl;
    cout << "commmand:" << endl;
    cout << "\tLCSCalculator -A MLCSAPP [-i input][-o output][-a alphabets][-e \"[k|c]\"]" << endl;
    cout << "someelse parameters:" << endl;
    cout << "\tk: The number of points maintained in the algorithm." << endl;
    cout << "\tc: The difference between the level of the best point and the level of the worst point." << endl;
    cout << endl;
}