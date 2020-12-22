#include "io.h"

MLCSIO::MLCSIO(string filename) : filename(filename){
	
	ifstream f(filename);
	if(f.good()){
		string seq;
		while(!f.eof()){
			getline(f, seq);
			if(seq.size() > 1) seqs.push_back(seq);
		}
	}
		
}

void MLCSIO::output(ostream& os, string algo, string alphasets){
	
	if(seqs.size() == 0){
		os << "Test Set is empty!\n";
	}
	else{
		os << "Test set:" << endl;
		for(int i = 0; i < seqs.size(); i++){
			os << seqs[i] << endl;
		}
	}
	if(algo == TOPMLCSSYM){
		TOP_MLCS topmlcs(seqs, alphasets);
		vector<string> mlcs;
		clock_t start_t, end_t;
		
		start_t = clock();
		topmlcs.run();
		end_t = clock();
		mlcs = topmlcs.MLCS();
		os << "Result(by " << algo << "):\n";
		os << "time(us) : " << end_t - start_t << "\n";
		os << "the length of lcs : " << mlcs[0].length() << "\n";
		os << "all lcs : \n";
		for(int i = 0; i < mlcs.size(); i++){
			os << i << " : " << mlcs[i] << "\n";
		}
	}
	else if(algo == MLCSAPPSYM){
		/* else argument*/
		int CONSTK, CONSTC;
		cout << "The const value k(negative value respresents default value) > ";
		cin >> CONSTK;
		cout << "The const value c(negative value respresents default value) > ";
		cin >> CONSTC;
		
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
	}
	else if(algo == PROMLCSSYM){
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
	}
	else if(algo == QUICKDPSYM){
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
	}
	else if(algo == RLPMLCSSYM){
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
	}
	else{
		os << "'" << algo << "' is undefined.\n";
	}
	
}
