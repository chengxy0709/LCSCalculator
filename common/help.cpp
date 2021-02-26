#include "help.h"

int UsageforAlg(string& alg){

    if(alg == HASMLCSSYM){
        UsageforHASMLCS();
    }
    else if(alg == MLCSAPPSYM){
        UsageforMLCSAPP();
    }
    else if(alg == PROMLCSSYM){
        UsageforPROMLCS();
    }
    else if(alg == TOPMLCSSYM){
        UsageforTOPMLCS();
    }
    else if(alg == QUICKDPSYM){
        UsageforQUICKDP();
    }
    else{
        return -1;
    }

    return 0;
}