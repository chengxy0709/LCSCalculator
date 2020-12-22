#include "MLCSAPP.h"

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

	ATTR(Attribute, p->attr)->f = h + ATTR(Attribute, p->attr)->g;
	
}


int MLCSAPP::getLCS( Point<CordType>* p ) {
	while(p->cord[0] != 0){
		lcs = seqs[0][p->cord[0] - 1] + lcs;
		p = ATTR(Attribute, p->attr)->parent;
	}
	return lcs.length();
}

struct Greater_Fun{
    size_t operator() (const Point<CordType>* a, const Point<CordType>* b) const{
        return ATTR(Attribute, a->attr)->f > ATTR(Attribute, b->attr)->f;
    }
};

void MLCSAPP::run(){

    HashSet nodeset(CONSTK);
	vector< vector< Point<CordType>* > > trash;
	Point<CordType> *p0 = new Point<CordType>(seqs.size(), false, 0);
	bool flag = true;
	
	p0->attr = COMADDR(new Attribute);
	computeF(p0);
	nodeset.insert(p0);

	while(nodeset.size() > 0) {
        vector< Point<CordType>* > nodeset_t;
        int maxf = 0;

        for(auto& p : nodeset){
            if(ATTR(Attribute, p->attr)->f > maxf) maxf = ATTR(Attribute, p->attr)->f;
        }
        for(auto& p : nodeset){
            if(ATTR(Attribute, p->attr)->f >= (maxf - CONSTC)){
                nodeset_t.push_back(p);
            }
			else{
				delete ATTR(Attribute, p->attr);
				delete p;
			}
        }

        if(nodeset_t.size() > CONSTK){
            sort(nodeset_t.begin(), nodeset_t.end(), Greater_Fun());
            for(unsigned int i = nodeset_t.size() - 1; i > CONSTK; i--){
                delete ATTR(Attribute, nodeset_t.back()->attr);
				delete nodeset_t.back();
                nodeset_t.pop_back();
            }
        }

		trash.push_back(nodeset_t);
        nodeset.clear();
        flag = true;
        for(auto& p : nodeset_t){
            if (ATTR(Attribute, p->attr)->f - ATTR(Attribute, p->attr)->g != 0) { // h(p) != 0
                flag = false;
            }

            for(int i = 0; i < SucTabs[0].size(); i++) {
				Point<CordType> *suc = successor(p, SucTabs, i);
				if(suc){
					auto q = nodeset.find(suc);
					if (q != nodeset.end()) { //node is already in node set
						Attribute *sucAttr = ATTR(Attribute, (*q)->attr);
						if ( sucAttr->g < ATTR(Attribute, p->attr)->g + 1) {
							sucAttr->parent = p;
							sucAttr->g = ATTR(Attribute, p->attr)->g + 1;
							computeF(*q);
						}
						delete suc;
					} 
					else {
						suc->attr = COMADDR(new Attribute(0, ATTR(Attribute, p->attr)->g + 1, p));
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
			delete ATTR(Attribute, p->attr);
			delete p;
		}
	}
	
}

MLCSAPP::~MLCSAPP(){
		
}