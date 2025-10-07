#ifndef LBUFFEREDMINHASH_HPP
#define LBUFFEREDMINHASH_HPP
#include <set>
#include <algorithm>
#include <iostream>
using namespace std;

class l_bufferd_minhash{
public:
    set<pair<int,int >> S;
    int l;
    pair<int,int >delta;
    void initialize(int num,pair<int,int> d){
        l=num,delta=d;
        S.clear();
        return ;
    }
    void insert(int hash_value,int u){
        if(pair<int,int>(hash_value,u)>delta) 
            return ;
        if((int)S.size()==l){
            if((*prev(S.end())).first>hash_value){
                S.insert(pair<int,int>(hash_value,u));
                S.erase(prev(S.end()));
            }
            delta=(*prev(S.end()));
        }
        else{
            S.insert(pair<int,int>(hash_value,u));
        }
        return ;
    }
    bool delete_element(int hash_value,int u){
        set<pair<int,int >>::iterator ele=S.find(pair<int,int>(hash_value,u));
        if(ele!=S.end()){
            S.erase(ele);
        }
        return S.empty();
    }
};

#endif