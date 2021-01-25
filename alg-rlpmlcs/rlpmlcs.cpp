#include "rlpmlcs.h"

static vector<int> firstlayer;

RLP_MLCS::RLP_MLCS(vector<string>& seqs, string& alphabets){

	unordered_map<char, int> cmap = build_alphabet_map(alphabets);
	phash_init(seqs.size());
	SucTabs = cal_suc_tabs(seqs, cmap);
	this->seqs = seqs;

}

int RLP_MLCS::cal_opt_point(){
    
    Point<CordType> *p0 = new Point<CordType>(SucTabs.size(), false, 0);
    int index = -1;
    int mincord = INT32_MAX;

    DM.insert(p0);
    for(int i = 0; i < SucTabs[0].size(); i++){
        Point<CordType> *suc = successor(p0, SucTabs, i);
        int cord = 0;
        if(suc){
            DM.insert(suc);
            firstlayer.push_back(DM.at(suc));
            for(int j = 0; j < SucTabs.size(); j++){
                cord += suc->cord[j];
            }
            if(cord < mincord){
                index = i;
                mincord = cord;
            }
        }
    }

    return index;
}

int RLP_MLCS::construct_ICSG(Point<CordType> *p0){

    queue< Point<CordType>* > Q;
    int index = 0;
	
    Q.push(p0);
    if(ID.size() < DM.size()) ID.resize(DM.size(), 0);

    while(!Q.empty()){
        Point<CordType> *p = Q.front();
        for(int i = 0; i < SucTabs[0].size(); ++i){
            Point<CordType> *suc = successor(p, SucTabs, i);
            if(suc){
                if(DM.insert(suc)){
                    Q.push(suc);
                    ID.push_back(1);
                }
                else{
                    int index_t = DM.at(suc);
                    ID[index_t] = ID[index_t] + 1;
					delete suc;
                }
            }
        }
        Q.pop();
        index++;
    }

    Point<CordType> *pinf = new Point<CordType>(SucTabs.size(), false, -1);
    DM.insert(pinf);

    return DM.size();

}

int RLP_MLCS::ForwardTopSort(int maxIndex, int startpoint, int curlevel){

    vector<int> D{startpoint};
    tlevel = vector<int>(maxIndex, 0);
    precursor = vector< list<int> >(maxIndex, list<int>());

    int k = curlevel;
    tlevel[startpoint] = curlevel;

    for(int i = 0; i < firstlayer.size(); i++){
        precursor[firstlayer[i]].push_back(0);
    }

    while(D.size()){
        vector<int> Dt;
        for(auto index : D){
            bool hasSuccessor = false;
			Point<CordType> *p = DM.at(index);
			for(int i = 0; i < SucTabs[0].size(); ++i){
				Point<CordType> *suc = successor(p, SucTabs, i);
                if(suc){
                    int index_t = DM.at(suc);
                    hasSuccessor = true;
                    tlevel[index_t] = k + 1;
                    if(is_immediate_successor(suc, p, SucTabs)){
                        precursor[index_t].push_back(index);
                    }
                    ID[index_t] = ID[index_t] - 1;
                    if(ID[index_t] == 0){
                        Dt.push_back(index_t);
                    }
					delete suc;
                }
            }
            if(!hasSuccessor){
                precursor[maxIndex - 1].push_back(index);
            }
        }
        k = k + 1;
        D = Dt;
    }
    tlevel[maxIndex - 1] = k;

    return k - 1;

}



void RLP_MLCS::BackwardTopSort(int maxlevel, int pinf_index){

    vector<int> D{pinf_index};
    int k = 0;

    inLCSPath = vector<bool>(pinf_index, false);
    while(D.size() != 0){
        vector<int> Dt;
        for(auto q : D){
            for(auto p = precursor[q].begin(); p != precursor[q].end();){
                if(tlevel[*p] + k != maxlevel){
                    p = precursor[q].erase(p);
                }
                else{
                    Dt.push_back(*p);
                    p++;
                }
            }
            if(precursor[q].size() > 0){
                inLCSPath[q] = true;
            }
        }
        D = Dt;
        k = k + 1;
    }

}

void RLP_MLCS::construct_SubICSG(int startpoint, int curlevel){

    queue<int> Q;

    tlevel[startpoint] = curlevel;
    Q.push(startpoint);
    if(tlevel.size() < DM.size()) tlevel.resize(DM.size());

    while(!Q.empty()){
        vector<int> Dt;
        Point<CordType> *p = DM.at(Q.front());
        for(int i = 0; i < SucTabs[0].size(); i++){
            Point<CordType> *suc = successor(p, SucTabs, i);
            if(suc){
                if(DM.insert(suc)){ // new point
                    tlevel.push_back(tlevel[Q.front()] + 1);
                    precursor.push_back({Q.front()});
                    Q.push(DM.size() - 1);
                }
                else{ // 
                    int index = DM.at(suc);
                    if(tlevel[index] < tlevel[Q.front()] + 1){
                        precursor[index].clear();
                        precursor[index].push_back(Q.front());
                        tlevel[index] = tlevel[Q.front()] + 1;
                    }
                    else if(tlevel[index] == tlevel[Q.front()] + 1){
                        precursor[index].push_back(Q.front());
                    }
                    delete suc;
                }
            }
        }
        Q.pop();
    }

}

void RLP_MLCS::run(){
	
	int maxIndex, maxlevel, optpind;
	string LCSRecord;
	
    optpind = cal_opt_point(); // optimized point's index in firstlayer
    maxIndex = construct_ICSG(DM.at(firstlayer[optpind]));
    cout << "The number of points Generated is " << maxIndex << " ." << endl;
    maxlevel = ForwardTopSort(maxIndex ,firstlayer[optpind] , 1);
    cout << "The maximal level of points is " << maxlevel << " ." << endl;
    BackwardTopSort(maxlevel, maxIndex - 1);
	cout << "Generating MLCS..." << endl;
    for(int i = 0; i < firstlayer.size(); i++){
        if(i == optpind) continue;
        construct_SubICSG(firstlayer[i], 1);
    }
    /*for(int i = 0; i < DM.size(); i++){
        cout << i << ':';
        DM.at(i)->print(SucTabs.size(), '\t');
        cout << tlevel[i] << '\t';
        cout << inLCSPath[i] << "\t|";
        for(auto ii : precursor[i]) cout << ii << ',';
        cout << endl;
    }*/
	
	GetMLCS(LCSRecord, maxIndex - 1);
	for(string& lcs : mlcs){
        reverse(lcs.begin(), lcs.end());
    }
	
}

void RLP_MLCS::GetMLCS(string& LCSRecord, int index){

    if(index == 0){
        mlcs.push_back(LCSRecord.substr(0, LCSRecord.length() - 1));
        return;
    }

    for(auto p : precursor[index]){
        Point<CordType> *q = DM.at(p);
        LCSRecord += seqs[0][q->cord[0] - 1];
        GetMLCS(LCSRecord, p);
        LCSRecord = LCSRecord.substr(0, LCSRecord.length() - 1);
    }

}

RLP_MLCS::~RLP_MLCS(){
	
	
}

int exe_rlpmlcs(vector<string>& seqs, string& alphasets, ostream& os, string& algo){
    if(algo == RLPMLCSSYM){
        RLP_MLCS rlpmlcs(seqs, alphasets);
		vector<string> mlcs;
		clock_t start_t, end_t;
		
		start_t = clock();
		rlpmlcs.run();
		end_t = clock();
		mlcs = rlpmlcs.MLCS();
		os << "Result(by " << algo << "):\n";
		os << "time(us) : " << end_t - start_t << "\n";
		os << "the length of lcs : " << mlcs[0].length() << "\n";
		os << "all lcs : \n";
		for(int i = 0; i < mlcs.size(); i++){
			os << i << " : " << mlcs[i] << "\n";
		}
        return 0;
    }
    else{
        return -1;
    }
}