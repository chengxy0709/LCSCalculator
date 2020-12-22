#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

#define OPTEXPR "p:l:n:N:a:"

const string PATH = "./example";
const int LENGTH = 50;
const int NUM = 1;
const int ALPHABETSNO = 1;
const string ALPHABETS[] = {"ACGT",   \
                            "ACDEFGHIKLMNPQRSTVWY"};

struct Arguments{
    string path;
    int len; 
    int num;
    int alphaNo;
    string alphabets;
};

void Usage(){
    cout << "Usage:" << endl;
    cout << "-p <path>:      -----  specify a file store sequences." << endl;
    cout << "-l <length>:    -----  the length of a sequence." << endl;
    cout << "-n <num>        -----  the number of sequences." << endl; 
    cout << "-N <No>         -----  specify a existing alphabet set[1-2]." << endl;
    cout << "-a <alphabets>  -----  specify a unique alphabet set, No must be 0." << endl;   
}

string genOneSeq(int len, int alphaNo, string& alphabets){
    string seq;
    string alphabets_t = (alphaNo == 0) ? alphabets : ALPHABETS[alphaNo - 1];
    int alphabetsz = alphabets_t.length();
    for(int i = 0; i < len; i++){
        seq += alphabets_t[rand() % alphabetsz];
    }

    return seq;
}

vector<string> genSeqs(int len, int num, int alphaNo, string& alphabets){
    vector<string> seqs;
    for(int i = 0; i < num; i++){
        seqs.push_back(genOneSeq(len, alphaNo, alphabets));
    }

    return seqs;
}

void getArgs(int argc, char* argv[], Arguments& args){
    
    int ch;
    args.alphabets = ALPHABETS[ALPHABETSNO];
    args.alphaNo = ALPHABETSNO;
    args.len = LENGTH;
    args.num = NUM;
    args.path = PATH;
    while((ch = getopt(argc, argv, OPTEXPR)) != -1){
        switch (ch)
        {
        case 'a':
            args.alphabets = optarg;
            break;
        case 'N':
            args.alphaNo = atoi(optarg);
            break;
        case 'n':
            args.num = atoi(optarg);
            break;
        case 'l':
            args.len = atoi(optarg);
            break;
        case 'p':
            args.path = optarg;
            break;
        default:
            Usage();
            break;
        }
    }

}


int main( int argc, char* argv[] ){

    Arguments args;
    vector<string> output;
    ofstream of;

    getArgs(argc, argv, args);
    srand(time(NULL));
    
    output = genSeqs(args.len, args.num, args.alphaNo, args.alphabets);
    of.open(args.path, ofstream::out | ofstream::trunc);
    
    if(of.is_open()){
        for(auto& seq : output){
            of << seq << endl;
        }    
    }
    else{
        cerr << "failed to open the output file." << endl;
    }

    return 0;

}