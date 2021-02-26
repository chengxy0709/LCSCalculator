#include <iostream>
#include <vector>
#include <malloc.h>
#include <unistd.h>
#include <time.h>
#include <fstream>

#include "io.h"
#include "help.h"

using namespace std;

struct Arguments
{
    string alg;
    string outputfile;
    string inputfile;
    string alphabets;
    string params;
};


static void Usage(){
    cout << "Usage:" << endl;
    cout << "-a <alphabets>    -----  specify a file store sequences." << endl;
    cout << "-A <algorithm>    -----  select a algorithm <TOPMLCS|MLCSAPP|PROMLCS|HASMLCS|QUICKDP>." << endl;
    cout << "-o <outputfile>   -----  specify a output file." << endl; 
    cout << "-i <inputfile>    -----  specify a input file." << endl;
    cout << "-e <extra param>  -----  specify some extra parameters for a specified algorithm." << endl;
    cout << "-h <help message> -----  print some help messages for a specified algorithm." << endl;
}

int getargs(int argc, char* argv[], Arguments& args);

int main(int argc, char* argv[]){
    
    Arguments args;
    MLCSIO mlcs;
    vector<string> seqs{"ACTAGTGC", "TGCTAGCA", "CATGCGAT"};

    cout << "**********************************\n";
    cout << "*                                *\n";
    cout << "*        LCS   Calculator        *\n";
    cout << "*                                *\n";
    cout << "**********************************\n";
    int ret = getargs(argc, argv, args);
    if(ret == -1){
        cout << "parameters error." << endl;
        return -1;
    }
    else if(ret != 0){
        return 0;
    }

    if(args.inputfile.length() == 0){
        mlcs = MLCSIO(seqs);
    }
    else
    {
        mlcs = MLCSIO(args.inputfile);
    }
    if(args.outputfile.length() == 0){
        mlcs.output(cout, args.alg, args.alphabets, args.params);
    }
    else{
        filebuf fb;
        fb.open(args.outputfile, ios::out | ios::trunc);
        ostream of(&fb);
        mlcs.output(of, args.alg, args.alphabets, args.params);
    }

    return 0;

}

int getargs(int argc, char* argv[], Arguments& args){
    
    int ch;
    int flag = 0;
    args.alphabets = "ACGT";
    args.alg = "TOPMLCS";
    args.inputfile = "";
    args.outputfile = "";
    args.params = "";

    while((ch = getopt(argc, argv, "a:A:i:o:e:h:")) != -1){
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
        case 'e':
            args.params = optarg;
            break;
        case 'h':{
            string s = optarg;
            if(UsageforAlg(s)){
                cout << "invalid algorithm.";
            }
            flag = -2;
            break;
        }
        default:
            Usage();
            flag = -1;
        }
    }

    return flag;

}