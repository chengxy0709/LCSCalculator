#include "hasmlcs.h"

// let k to be 2 (It is 0.5 in article, but it is doubted for me to this value.)
#define KNORM 2
double HASMLCS::cal_k_norm(const Point<CordType>* p){
    double res = 0;
    for(int i = 0; i < SucTabs.size(); i++){
        res += pow(seqs[i].length() - p->cord[i], KNORM);
    }
    return pow(res, 1.0 / KNORM);
}

// like lexicographical order
bool porder(const Point<CordType>* p1, const Point<CordType>* p2){

    for(int i = 0; i < g_point_size; i++){
        if(p1->cord[i] < p2->cord[i]) return true;
        else if(p1->cord[i] > p2->cord[i]) return false;
    }

    // it is impossible for two different points to come here.
    return false;
    
}

bool priority1::comp(const Point<CordType>* p1, const Point<CordType>* p2) const{
    if(ATTR(HASAttr, p1)->f < ATTR(HASAttr, p2)->f){
        return true;
    }
    // same priority, do these steps to break this tie
    else if(ATTR(HASAttr, p1)->f == ATTR(HASAttr, p2)->f){
        // using level infomation
        if(ATTR(HASAttr, p1)->g < ATTR(HASAttr, p2)->g){
            return true; // bigger level
        }
        else if(ATTR(HASAttr, p1)->g > ATTR(HASAttr, p2)->g){
             return false;
        }
        else{ // using the second method, k(v)
            //return ATTR(HASAttr, p1)->k < ATTR(HASAttr, p2)->k;
            if(ATTR(HASAttr, p1)->k < ATTR(HASAttr, p2)->k) return true;
            else if(ATTR(HASAttr, p1)->k > ATTR(HASAttr, p2)->k) return false; 
            else return porder(p1, p2);
        }
    }
    else{
        return false;
    }
}

bool priority3::comp(const Point<CordType>* p1, const Point<CordType>* p2) const{
    if(ATTR(HASAttr, p1)->f > ATTR(HASAttr, p2)->f){
        return true;
    }
    // same priority, do these steps to break this tie
    else if(ATTR(HASAttr, p1)->f == ATTR(HASAttr, p2)->f){
        // using level infomation
        if(ATTR(HASAttr, p1)->g > ATTR(HASAttr, p2)->g){
            return true; // bigger level
        }
        else if(ATTR(HASAttr, p1)->g < ATTR(HASAttr, p2)->g){
             return false;
        }
        else{ // using the second method, k(v)
            //return ATTR(HASAttr, p1)->k < ATTR(HASAttr, p2)->k;
            if(ATTR(HASAttr, p1)->k > ATTR(HASAttr, p2)->k) return true;
            else if(ATTR(HASAttr, p1)->k < ATTR(HASAttr, p2)->k) return false; 
            else return porder(p1, p2);
        }
    }
    else{
        return false;
    }
}

const double Po[1][1] = {0};
double HASMLCS::cal_ex_val(const Point<CordType>* p){
    double lmin = INT32_MAX;
    double res = 0;
    for(int i = 0; i < SucTabs.size(); i++){
        if(lmin > seqs[i].length() - p->cord[i] + 1)
            lmin = seqs[i].length() - p->cord[i] + 1;
    }
    for(int i = 1; i <= lmin; i++){
        double po = 1;
        for(int j = 0; j < seqs.size(); j++){
            po *= Po[i - 1][seqs[j].length() - p->cord[j] + 1];
        }
        res += pow(1 - po, pow(seqs.size(), i));
    }
    return lmin - res;
}

bool priority2::comp(const Point<CordType>* p1, const Point<CordType>* p2) const{
    if(ATTR(HASAttr, p1)->EX > ATTR(HASAttr, p2)->EX){
        return true;
    }
    else{
        return false;
    }
}
