#ifndef IO_H
#define IO_H

#include <ostream>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

#include "TOP_MLCS.h"
#include "MLCSAPP.h"
#include "PRO_MLCS.h"
#include "QuickDP.h"

#define TOPMLCSSYM "TOPMLCS"
#define MLCSAPPSYM "MLCSAPP"
#define PROMLCSSYM "PROMLCS"
#define QUICKDPSYM "QUICKDP"

using namespace std;

class MLCSIO{
	
public:
	MLCSIO(){}
	MLCSIO(vector<string>& seqs):seqs(seqs){}
	MLCSIO(string filename);
	~MLCSIO(){}
	
	vector<string> getSeqs() const {return seqs;}
	void output(ostream& os, string algo, string alphasets);
	
private:
	vector<string> seqs;
	string filename;

};

#endif // IO_H