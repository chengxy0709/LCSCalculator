#include <iostream>
#include <vector>
#include <malloc.h>
#include <time.h>
#include <fstream>

#include "io.h"

using namespace std;

// argument list:
// algorithm
// algorithm outposition
// algorithm filename alphabets
// algorithm outposition filename alphabets
int main(int argc, char* argv[]){
	
	string str = "ATCG";
	vector<string> seqs{"ACTAGTGC", "TGCTAGCA", "CATGCGAT"};

	if(argc <= 1){
		cout << "argument too few.\n";
	}
	else if(argc == 2){
		MLCSIO mlcs(seqs);
		mlcs.output(cout, string(argv[1]), str);
	}
	else if(argc == 3){
		MLCSIO mlcs(seqs);
		filebuf fb;
		fb.open(string(argv[2]), ios::out);
		ostream of(&fb);
		mlcs.output(of, string(argv[1]), str);
	}
	else if(argc == 4){
		string filename(argv[2]);
		string alphabets(argv[3]);
		MLCSIO mlcs(filename);
		mlcs.output(cout, string(argv[1]), alphabets);
	}
	else if(argc == 5){
		string filename(argv[3]);
		string alphabets(argv[4]);
		MLCSIO mlcs(filename);
		filebuf fb;
		fb.open(string(argv[2]), ios::out);
		ostream of(&fb);
		mlcs.output(of, string(argv[1]), alphabets);
	}
	else{
		cout << "argument too many.\n";
	}

	return 0;

}
