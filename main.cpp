#include <iostream>
#include <vector>
#include <malloc.h>
#include <unistd.h>
#include <time.h>
#include <fstream>

#include "io.h"

using namespace std;

struct Arguments
{
	string alg;
	string outputfile;
	string inputfile;
	string alphabets;
};


static void Usage(){
    cout << "Usage:" << endl;
    cout << "-a <alphabets>    -----  specify a file store sequences." << endl;
    cout << "-A <algorithm>    -----  select a algorithm <TOPMLCS|MLCSAPP|PROMLCS|QUICKDP>." << endl;
    cout << "-o <outputfile>   -----  specify a output file." << endl; 
    cout << "-i <inputfile>    -----  specify a input file." << endl;
	cout << "\nusing some default options..." << endl;
}

void getargs(int argc, char* argv[], Arguments& args);

int main(int argc, char* argv[]){
	
	Arguments args;
	MLCSIO mlcs;
	vector<string> seqs{"ACTAGTGC", "TGCTAGCA", "CATGCGAT"};

	cout << "LCS Calculator\n";
	getargs(argc, argv, args);
	if(args.inputfile.length() == 0){
		mlcs = MLCSIO(seqs);
	}
	else
	{
		mlcs = MLCSIO(args.inputfile);
	}
	if(args.outputfile.length() == 0){
		mlcs.output(cout, args.alg, args.alphabets);
	}
	else{
		filebuf fb;
		fb.open(args.outputfile, ios::out | ios::trunc);
		ostream of(&fb);
		mlcs.output(of, args.alg, args.alphabets);
	}

	return 0;

}

void getargs(int argc, char* argv[], Arguments& args){
    
    int ch;
    args.alphabets = "ACGT";
    args.alg = "TOPMLCS";
    args.inputfile = "";
    args.outputfile = "";

    while((ch = getopt(argc, argv, "a:A:i:o:")) != -1){
        switch (ch)
        {
        case 'a':
            args.alphabets = optarg;
            break;
        case 'A':
            args.alg = optarg;
            break;
        case 'i':
            args.inputfile = optarg;
            break;
        case 'o':
            args.outputfile = optarg;
            break;
        default:
            Usage();
            break;
        }
    }

}