#include "promlcs.h"

struct comp{
    bool operator() (const Point<CordType> *p1, const Point<CordType> *p2){
        return ATTRINT(p1) > ATTRINT(p2);
    }
};

class Layer{

public:
    Layer(int d){dtree.setD(d);order = 0;}
    ~Layer(){
        while(!priq.empty()){
            delete priq.top();
            priq.pop();
        }
    }

    void setOrder(int ord){order = ord;}
    int getOrder() const {return order;}
    void addpriq(Point<CordType> *p){priq.push(p);}
    Point<CordType>* rmpriq(){
        Point<CordType> *p = priq.top();
        priq.pop();
        return p;
    }
    bool addDtree(Point<CordType> *p){
        if(dtree.Search(p, dtree.getRoot())){
            dtree.Insert(p, dtree.getRoot());
            return true;
        }
        else return false;
    }
    bool isPriqEmpty(){return priq.empty();}
    void printInfo(){
        cout << priq.size() << endl;
        //dtree.traverse(dtree.getRoot());
    }

private:
    priority_queue<Point<CordType>*, vector< Point<CordType>* >, comp> priq;
    int order;
    Dtree dtree;

};

PRO_MLCS::PRO_MLCS(vector<string> seqs, string alphabets, int scnum)
 : seqs(seqs),alphaSize(alphabets.length())
{

    unordered_map<char, int> cmap = build_alphabet_map(alphabets);
    if(scnum > 0) this->scnum = scnum;
    else scnum = 100;
    SucTabs = cal_suc_tabs(seqs, cmap);
    calSU(cmap);
    
}

void PRO_MLCS::calSU(unordered_map<char, int>& cmap){

    for(int i = 0; i < seqs.size(); i++){
        int len = seqs[i].length();
        SU.push_back(vector< vector<int> >(alphaSize, vector<int>(len + 1, 0)));
        // calculate SU[i]
        for (int j = len - 1; j >= 0; j--) {
            for (int k = 0; k < alphaSize; k++) {
                SU[i][k][j] = SU[i][k][j + 1];
            }
            SU[i][cmap[seqs[i][j]]][j]++;
        }
    }
    
}

int PRO_MLCS::calg(Point<CordType> *p, int order){

    int f = 0;
    int minv = 0;
    for(int i = 0; i < alphaSize; i++){
        minv = INT_MAX;
        for(int j = 0; j < seqs.size(); j++){
            if(minv > SU[j][i][p->cord[j]]){
                minv = SU[j][i][p->cord[j]];
            }
        }
        f = minv + f;
    }

    return f + order;

}

int PRO_MLCS::caldist(Point<CordType> *p){
    int sum = 0;
    for(int i = 0; i < seqs.size(); ++i){
        sum += p->cord[i];
    }
    return sum;
}

void PRO_MLCS::run(){
    
    list<Layer*> layers;
    list<Layer*>::iterator lbeginIter;
    list<Layer*>::iterator lendIter;
    list<Layer*>::iterator curIter;

    int Lapp = 0;

    Point<CordType> *p = new Point<CordType>(seqs.size(), false, 0);
    Layer *layer = new Layer(seqs.size());

    SETATTRINT(p, 0);
    layer->addpriq(p);

    layers.push_back(layer);
    lbeginIter = layers.begin();
    curIter = layers.begin();
    lendIter = layers.begin();

    Lapp = 0;

    while(1){
        int m = 0;
        Layer *cur = *curIter;
        auto cur_next = curIter;
        cur_next++;
        while(!cur->isPriqEmpty()){
            p = cur->rmpriq();
            if(calg(p, cur->getOrder()) > Lapp && cur->addDtree(p)){
                m++;
                for(int i = 0; i < alphaSize; ++i){
                    Point<CordType> *suc = successor(p, SucTabs, i);
                    if(suc){
                        if(curIter == lendIter){
                            layer = new Layer(seqs.size());
                            layer->setOrder(cur->getOrder() + 1);
                            layers.push_back(layer);
                            Lapp++;
                            lendIter++;
                            cur_next = lendIter;
                        }
                        (*cur_next)->addpriq(suc);
                    }
                }
            }
            else{
                delete p;    // be careful of this statement.
                            // Don't delete a point that has already been stored in d-tree.
            }
            if(m == scnum){
                break;
            }
        }
        if(lbeginIter == lendIter){
            break;
        }
        if(curIter == lbeginIter && cur->isPriqEmpty()){
            lbeginIter++;
        }
        if(curIter != lendIter){
            curIter++;
        }
        else{
            break;
            curIter = lbeginIter;
        }
//        cout << "cur: " << (*curIter)->getOrder() << endl;
//        cout << "begin: " << (*lbeginIter)->getOrder() << endl;
//        cout << "end: " << (*lendIter)->getOrder() << endl;
    }

    cout << "The length of LCS is " << Lapp << " ." << endl;
    
    // free memory
    for(auto& l : layers){
        delete l;
    }
    
}

int exe_promlcs(vector<string>& seqs, string& alphasets, ostream& os, string& algo){
    if(algo == PROMLCSSYM){
        /* else argument*/
        int scnum;
        cout << "The number of calculating point in a single iteration(negative value respresents default value) > ";
        cin >> scnum;
        
        PRO_MLCS promlcs(seqs, alphasets, scnum);
        string lcs;
        clock_t start_t, end_t;
        
        start_t = clock();
        promlcs.run();
        end_t = clock();
        lcs = promlcs.LCS();
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