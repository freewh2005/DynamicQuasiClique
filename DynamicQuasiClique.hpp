#ifndef DYNAMICQUASICLIQUE_HPP
#define DYNAMICQUASICLIQUE_HPP
#include "getArgs.hpp"
#include "LBufferedMinHash.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <string>
#include <unordered_map>
using namespace std;
int n,m,use_bottom_k;
int operations = 0;
int k,l = 1;
int p = pow(2,31) - 1;
double ts;
double f;
double b;
double alpha;
int q;
int B = 5;
int Batch = 50000;
int QLimited = 100000;
const double sol_similarity = 0.7;
const double sol_admitted = 0.8;
ifstream fin;

vector<vector<int>> neighbor;
vector<int> deg;
vector<int> deg_gamma;
vector<pair<int,int>> sorted_vertex;

struct pair_hash {
    size_t operator()(const pair<int, int>& p) const {
        size_t seed = hash<int>()(p.first);
        seed ^= hash<int>()(p.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

unordered_map<pair<int,int>,bool,pair_hash>mp;
int sol_cnt;
int o,e,cnt;
vector<int> vec_a;
vector<int> vec_b;
vector<int> vis;
vector<vector<l_bufferd_minhash>> hash_set;

vector<set<pair<unsigned int ,int >> >Sk;
vector<vector<pair<unsigned int ,int >> >bottom_k;
vector<int >lazy;
vector<unsigned int>hash_bottom;

int deg_gamma_min;
int qccnt;
unsigned long long rand_a,rand_b;
unsigned long long randull(){
    unsigned long long r1 = (unsigned long long)(rand() & 0x7FFF);
    unsigned long long r2 = (unsigned long long)(rand() & 0x7FFF);
    unsigned long long r3 = (unsigned long long)(rand() & 0x7FFF);
    unsigned long long r4 = (unsigned long long)(rand() & 0x7FFF);
    unsigned long long r5 = (unsigned long long)(rand() & 0xF);          
    return (r1 << 49) | (r2 << 34) | (r3 << 19) | (r4 << 4) | r5;
}
void initial_variable(){
    deg.resize(n);
    deg_gamma.resize(n);
    neighbor.resize(n);
    vis.resize(n);
    if(!use_bottom_k){
        vec_a.resize(k);
        vec_b.resize(k);
        hash_set.resize(n);
        for(int i=0;i<k;i++){
            vec_a[i]=rand()%(p-1)+1;
            vec_b[i]=rand();
        }
        for(int i=0;i<n;i++){
            hash_set[i].resize(k);
        }
    }
    else{
        rand_a=randull();
        rand_b=randull();
        Sk.resize(n);
        bottom_k.resize(n);
        lazy.resize(n);
        hash_bottom.resize(n);
        for(int i=0;i<n;i++){
            bottom_k[i].resize(k);
            lazy[i]=0;
            hash_bottom[i]=(rand_a*i+rand_b)>>32;
        }
    }
    return ;
}
void read_graph(argsController &aC){
    use_bottom_k=0;
    string filePath;
    if (aC.exist("-f"))
        filePath = aC.get("-f");
    else {
        cout << "-f filePath" << endl;
        exit(-1);
    }
    if (aC.exist("-gamma"))
        ts = stod(aC.get("-gamma").c_str());

    if (aC.exist("-b"))
        b = stod(aC.get("-b").c_str());

    if (aC.exist("-k"))
        k = atoi(aC.get("-k").c_str());

    if (aC.exist("-l"))
        l = atoi(aC.get("-l").c_str());
    
    if (aC.exist("-alpha"))
        alpha = stod(aC.get("-alpha").c_str());

    if (aC.exist("-use_bottom_k"))
        use_bottom_k=1;

    if (aC.exist("-B"))
        B = atoi(aC.get("-B").c_str());

    if (aC.exist("-Batch"))
        Batch = atoi(aC.get("-Batch").c_str());
    
    if (aC.exist("-Q"))
        QLimited = atoi(aC.get("-Q").c_str());

    fin.open(filePath);
    fin>>n>>m>>q;
    int u,v;
    initial_variable();
    for(int i=0;i<m;i++){
        fin>>u>>v;
        deg[u]++,deg[v]++;
        neighbor[u].emplace_back(v);
        neighbor[v].emplace_back(u);
        mp[make_pair(u,v)]=1;
        mp[make_pair(v,u)]=1;
    }
}
int core_num(int u){
    int num=0;
    for(int v:neighbor[u]){
        if((deg[v]+1)>=(deg[u]+1)*ts) num++;
    }
    return num;
}
void sort_vertex(){
    for(int i=0;i<n;i++){
        deg_gamma[i]=core_num(i);
        if(deg_gamma[i]>0) sorted_vertex.emplace_back(pair<int,int>(deg_gamma[i],i));
    }
    sort(sorted_vertex.begin(),sorted_vertex.end(),greater<pair<int,int>>());
    return ;
}
void get_k_signatures(int u){
    if(!use_bottom_k){
        hash_set[u].clear();
        for(int i=0;i<k;i++){
            l_bufferd_minhash l_buffered;
            l_buffered.initialize(l,pair<int,int>(p,p));
            for(int j:neighbor[u]){
                l_buffered.insert(int(((long long)vec_a[i]*(long long)(j)+(long long)vec_b[i])%p),j);
            }
            l_buffered.insert(int(((long long)vec_a[i]*(long long)(u)+(long long)vec_b[i])%p),u);
            hash_set[u].emplace_back(l_buffered);
        }
    }
    else{
        Sk[u].clear();
        for(int j:neighbor[u]){
            Sk[u].insert(make_pair(hash_bottom[j],j));
        }
        Sk[u].insert(make_pair(hash_bottom[u],u));
        int siz=min(k,(int)Sk[u].size());
        set<pair<unsigned int ,int >>::iterator it=Sk[u].begin();
        for(int i=0;i<siz;i++){
            bottom_k[u][i]=(*it);
            it++;
        }
        lazy[u]=0;
    }
    return ;
}
void remakebottom_k(int u){
    if(!vis[u]){
        get_k_signatures(u);
        vis[u]=1;
    }
    int siz=min(k,(int)Sk[u].size());
    set<pair<unsigned int ,int >>::iterator it=Sk[u].begin();
    for(int i=0;i<siz;i++){
        bottom_k[u][i]=(*it);
        it++;
    }
    lazy[u]=0;
    return ;
}
vector<pair<unsigned int ,int > >tmp_vec;
double ct_score_minhash(int u,int v,int jarc=0){
    if(!use_bottom_k){
        int res=0;
        for(int i=0;i<k;i++) 
            if((*hash_set[u][i].S.begin()).first==(*hash_set[v][i].S.begin()).first) res++;
        double j_sim=(double)res/(double)k;
        if(jarc) return j_sim;
        return (double)(deg[v]+deg[u]+2)*j_sim/(j_sim+1)/(double)(deg[u]+1);
    }
    else{
        if(lazy[u]||!vis[u]) remakebottom_k(u);
        if(lazy[v]||!vis[v]) remakebottom_k(v);
        int sizu=min(k,(int)Sk[u].size()),sizv=min(k,(int)Sk[v].size());
        tmp_vec.clear();
        tmp_vec.resize(sizu+sizv);
        for(int i=0;i<sizu;i++)
            tmp_vec[i]=bottom_k[u][i];
        for(int i=0;i<sizv;i++)
            tmp_vec[sizu+i]=bottom_k[v][i];
        sort(tmp_vec.begin(),tmp_vec.end());
        int num=0,res=0;
        for(int i=0;i<sizu+sizv;){
            num++;
            if(i+1<sizu+sizv&&tmp_vec[i].second==tmp_vec[i+1].second) res++,i+=2;
            else i+=1;
            if(num==k) break;
        }
        double j_sim=(double)res/(double)k;
        if(jarc) return j_sim;
        return (double)(deg[v]+deg[u]+2)*j_sim/(j_sim+1)/(double)(deg[u]+1);
    }
}
int sift_num(int u){
    qccnt++;
    int res=0;
    if(vis[u]==0){
        get_k_signatures(u);
        vis[u]=1;
    }
    for(int v:neighbor[u]){
        if(!vis[v]){
            get_k_signatures(v);
            vis[v]=1;
        }
        if(ct_score_minhash(u,v)>=ts) res++;
    }
    double tmp=(double)res/(double)(deg[u]+1);
    if(tmp<b) return 0;
    return res+1;
}
void in_case_empty(int u,int i,int x,int v){
    if(hash_set[u][i].delete_element(x,v)){
        hash_set[u][i].initialize(l,pair<int,int>(p,p));
        for(int j:neighbor[u]){
            hash_set[u][i].insert(int(((long long)vec_a[i]*(long long)(j)+(long long)vec_b[i])%p),j);
        }
        hash_set[u][i].insert(int(((long long)vec_a[i]*(long long)(u)+(long long)vec_b[i])%p),u);
    }
    return ;
}
void remove_from_set(int u,unsigned int x,int v){
    if(!vis[u]) return ;
    auto it=Sk[u].find(make_pair(x,v));
    if(it!=Sk[u].end()){
        lazy[u]=lazy[u]|((int)Sk[u].size()<=k||x<=bottom_k[u][k-1].first);
        Sk[u].erase(it);
        return ;
    }
}

struct solution{
    int cur_sol_vertex = 0;
    int cur_sol_size = 0;
    int cur_sol_edges = 0;
    double cur_sol_quasi = 0;
    vector<int> cur_sol;
    vector<int> cur_deg;
    vector<int> cur_vis;
    int flag;
    void initial(){
        cur_deg.resize(n,0);
        cur_vis.resize(n,0);
        cur_sol_size=0;
        cur_sol_edges=0;
        flag=0;
    }
    void make_cur_sol(){
        cur_sol.clear();
        cur_sol_edges=0;
        cur_sol_size=1;
        cur_sol.emplace_back(cur_sol_vertex);
        for(int v:neighbor[cur_sol_vertex]){
            cur_deg[v]=0;
            if(ct_score_minhash(cur_sol_vertex,v)>ts) cur_sol.emplace_back(v),cur_sol_size++;
        }
        return ;
    }
    void make_cur_quasi(){
        flag++;
        for(int v:cur_sol) cur_vis[v]=flag;
        int unvis_cnt=0;
        for(int i:cur_sol){
            int unvis_for_i=cur_sol_size-1;
            for(int j:cur_sol){
                if(j!=i&&mp[make_pair(i,j)]) cur_sol_edges++,cur_deg[i]++,unvis_for_i--;
            }
            unvis_cnt+=unvis_for_i;
            if(unvis_cnt>(1-sol_admitted)*(double)cur_sol_size*(double)(cur_sol_size-1)){
                cur_sol_quasi=0;
                return ;
            }
        }
        if(cur_sol_size>1) cur_sol_quasi=double(cur_sol_edges)/double(cur_sol_size*(cur_sol_size-1));
        else cur_sol_quasi=1;
        return ;
    }
    void update_sol(){
        flag++;
        for(int v:cur_sol){
            cur_deg[v]=0;
        }
        cur_sol_edges=0;
        for(int v:neighbor[cur_sol_vertex]) cur_deg[v]=0;
        for(int v:cur_sol) cur_vis[v]=flag;
        int unvis_cnt=0;
        for(int i:cur_sol){
            int unvis_for_i=cur_sol_size-1;
            for(int j:cur_sol){
                if(j!=i&&mp[make_pair(i,j)]) cur_sol_edges++,cur_deg[i]++,unvis_for_i--;
            }
            unvis_cnt+=unvis_for_i;
            if(unvis_cnt>0.2*(double)cur_sol_size*(double)(cur_sol_size-1)){
                cur_sol_quasi=0;
                return ;
            }
        }
        if(cur_sol_size>1) cur_sol_quasi=double(cur_sol_edges)/double(cur_sol_size*(cur_sol_size-1));
        else cur_sol_quasi=1;
        return ;
    }
    bool WarningForQuasiNoEquated(){
        int edges=0;
        for(int i:cur_sol){
            for(int j:neighbor[i]){
                if(find(cur_sol.begin(),cur_sol.end(),j)!=cur_sol.end()){
                    edges++;
                }
            }
        }
        int siz=(int)cur_sol.size();
        double sol_quasi;
        if(siz>1) sol_quasi=double(edges)/double(siz*(siz-1));
        else sol_quasi=1;
        return sol_quasi!=cur_sol_quasi;
    }
    void add(int u,int v){
        if(cur_deg[u]&&cur_deg[v]){
            cur_deg[u]++,cur_deg[v]++;
            cur_sol_edges+=2;
            cur_sol_quasi=double(cur_sol_edges)/double(cur_sol_size*(cur_sol_size-1));
        }
        return ;
    }
    bool choose_one(int u){
        double next_sol_quasi=double(cur_sol_edges-cur_deg[u]*2)/double((cur_sol_size-2)*(cur_sol_size-1));
        if(next_sol_quasi>=alpha){
            cur_sol_size--;
            cur_sol_edges-=cur_deg[u]*2;
            cur_deg[u]=0;
            for(int v:neighbor[u])
                if(cur_deg[v]>0) cur_deg[v]--;
            cur_sol.erase(remove(cur_sol.begin(),cur_sol.end(),u),cur_sol.end());
            cur_sol_quasi=next_sol_quasi;
            o++;
            return 1;
        }
        return 0;
    }
    bool delete_one_element(int u,int v,int w){
        if(choose_one(u)) return 1;
        if(choose_one(v)) return 1;
        if(choose_one(w)) return 1;
        return 0;
    }

    bool choose_two(int u,int v){
        if(u==-1||v==-1) return 0;
        if(u==v) return 0;
        if(cur_sol_size==2) return 0;
        auto find_u_in_v_neighbor=find(neighbor[u].begin(),neighbor[u].end(),v);
        int inc=2*(find_u_in_v_neighbor!=neighbor[u].end());
        double next_sol_quasi=double(cur_sol_edges+inc-cur_deg[u]*2-cur_deg[v]*2)/double((cur_sol_size-3)*(cur_sol_size-2));
        if(next_sol_quasi>=alpha){
            cur_sol_size-=2;
            cur_sol_edges-=(cur_deg[u]*2+cur_deg[v]*2-inc);
            for(int w:neighbor[u])
                if(cur_deg[w]>0) cur_deg[w]--;
            for(int w:neighbor[v])
                if(cur_deg[w]>0) cur_deg[w]--;
            cur_deg[u]=0,cur_deg[v]=0;
            cur_sol.erase(remove(cur_sol.begin(),cur_sol.end(),u),cur_sol.end());
            cur_sol.erase(remove(cur_sol.begin(),cur_sol.end(),v),cur_sol.end());
            cur_sol_quasi=next_sol_quasi;
            e++;
            return 1;
        }
        return 0;
    }

    bool delete_two_element(int u,int v,int w,int z){
        if(choose_two(u,v)) return 1;
        if(choose_two(u,w)) return 1;
        if(choose_two(u,z)) return 1;
        if(choose_two(v,w)) return 1;
        if(choose_two(v,z)) return 1;
        if(choose_two(w,z)) return 1;
        return 0;
    }
    bool delete_zero_element(int u,int v){
        auto find_u_in_sol=find(cur_sol.begin(),cur_sol.end(),u);
        auto find_v_in_sol=find(cur_sol.begin(),cur_sol.end(),v);
        if(find_u_in_sol==cur_sol.end()||find_v_in_sol==cur_sol.end()) return 1;
        cur_deg[u]--,cur_deg[v]--;
        cur_sol_edges-=2;
        if(cur_deg[u]&&cur_deg[v]){
            cur_sol_quasi=double(cur_sol_edges)/double(cur_sol_size*(cur_sol_size-1));
            if(cur_sol_quasi>=alpha) return 1;
        }
        return 0;
    }
    int delete_how_much(int u,int v){
        if(delete_zero_element(u,v)) return 0;
        int A=-1,B=-1;
        for(int x:cur_sol){
            if(A==-1||deg[x]<deg[A]){B=A;A=x;}
            else if(B==-1||deg[x]<deg[B]) B=x;
        }
        if(delete_one_element(u,v,A)) return 1;
        if(delete_two_element(u,v,A,B)) return 2;
        return 3;
    }
};

vector<int >sol_vis;
int cp_flag = 0;
double solution_compare(const solution &u, const solution &v) {
    cp_flag++;
    for (int x : u.cur_sol) {
        sol_vis[x] = cp_flag;
    }
    int cnt = 0;
    for (int x : v.cur_sol) {
        if (sol_vis[x] == cp_flag) cnt++;
    }
    return (double)cnt / (double)(u.cur_sol_size + v.cur_sol_size - cnt);
}

vector<solution > sol;
int insert_sol(int u,int num){
    if(!num&&sol_cnt==B) return 0;
    int replace=-1,replacemin=-1;
    for(int i=0;i<sol_cnt;i++)
        if(u==sol[i].cur_sol_vertex) return 0;
    solution tmp;
    tmp.cur_sol_vertex=u;
    tmp.cur_sol_size=num;
    tmp.initial();
    tmp.make_cur_sol();
    if(sol_cnt<B){
        for(int i=0;i<sol_cnt;i++){
            if(solution_compare(sol[i],tmp)>sol_similarity) return 0;
        }
        tmp.make_cur_quasi();
        if(tmp.cur_sol_quasi<sol_admitted) return 0;
        sol[sol_cnt]=tmp;
        sol_cnt++;
        return 1;
    }
    for(int i=0;i<sol_cnt;i++){
        if(num>sol[i].cur_sol_size&&solution_compare(sol[i],tmp)>sol_similarity){
            replace=(replace==-1||sol[i].cur_sol_size<sol[replace].cur_sol_size)?i:replace;
        }
        if(num>sol[i].cur_sol_size){
            replacemin=(replacemin==-1||sol[i].cur_sol_size<sol[replacemin].cur_sol_size)?i:replacemin;
        }
    }
    if(replacemin==-1) return 0;
    tmp.make_cur_quasi();
    if(tmp.cur_sol_quasi<sol_admitted) return 0;
    if(replace!=-1){
        if(tmp.cur_sol_quasi>=min(alpha,sol[replace].cur_sol_quasi)){
            sol[replace]=tmp;
            return 1;
        }
        return 0;
    }
    if(tmp.cur_sol_quasi>=min(alpha,sol[replacemin].cur_sol_quasi)){
        sol[replacemin]=tmp;
        return 1;
    }
    return 0;
}
void iterate_vertex(){
    for(int i=0;i<(int)sorted_vertex.size();i++){
        auto u=sorted_vertex[i];
        int flag=(sol_cnt==0);
        for(int j=0;j<sol_cnt;j++){
            flag|=(u.first>=sol[j].cur_sol_size);
        }
        if(!flag&&sol_cnt==B) break;
        int num=sift_num(u.second);
        insert_sol(u.second,num);
    }
    return ;
}
void addNewNodeIntoKMinHash(int u,int v){
    if(!vis[u]){
        get_k_signatures(u);
        vis[u]=1;
    }
    for(int i=0;i<k;i++){
        hash_set[u][i].insert(int(((long long)vec_a[i]*(long long)(v)+(long long)vec_b[i])%p),v);
    }
    return ;
}
void addNewNodeIntoBottomK(int u,int v){
    if(!vis[u]){
        get_k_signatures(u);
        vis[u]=1;
    }
    lazy[u]=lazy[u]|((int)Sk[u].size()<k||bottom_k[u][k-1].first>hash_bottom[v]);
    Sk[u].insert(make_pair(hash_bottom[v],v));
    return ;
}
void add_edge(int u,int v){
    mp[make_pair(u,v)]=1;
    mp[make_pair(v,u)]=1;
    deg[u]++,deg[v]++;
    if(!use_bottom_k){
        addNewNodeIntoKMinHash(u,v);
        addNewNodeIntoKMinHash(v,u);
    }
    else{
        addNewNodeIntoBottomK(u,v);
        addNewNodeIntoBottomK(v,u);
    }
    neighbor[u].emplace_back(v);
    neighbor[v].emplace_back(u);
    int flag=0;
    for(int w:neighbor[u]){
        if((deg[u]+1)>=(deg[w]+1)*ts&&(deg[u]<(deg[w]+1)*ts||w==v)){
            deg_gamma[w]++;
            if(deg_gamma[w]>=deg_gamma_min){
                int num=sift_num(w);
                flag|=insert_sol(w,num);
            }
        }
    }
    for(int w:neighbor[v]){
        if((deg[v]+1)>=(deg[w]+1)*ts&&(deg[v]<(deg[w]+1)*ts||w==u)){
            deg_gamma[w]++;
            if(deg_gamma[w]>=deg_gamma_min){
                int num=sift_num(w);
                flag|=insert_sol(w,num);
            }
        }
    }
    if(flag){
        for(int i=0;i<sol_cnt;i++)
            sol[i].update_sol();
    }
    else{
        for(int i=0;i<sol_cnt;i++){
            sol[i].add(u,v);
        }
    }
    return ;
}
void initial_request_for_ans(){
    int ans=-1;
    deg_gamma_min=n;
    for(int i=0;i<sol_cnt;i++){
        sol[i].update_sol();
        if(sol[i].cur_sol_quasi>=alpha){
            if(ans==-1||sol[i].cur_sol_size>sol[ans].cur_sol_size){
                ans=i;
            }
        }
        deg_gamma_min=min(deg_gamma_min,sol[i].cur_sol_size);
    }
    if(ans==-1){
        for(int i=0;i<sol_cnt;i++){
            if(ans==-1||sol[i].cur_sol_quasi>=sol[ans].cur_sol_quasi){
                ans=i;
            }
        }
    }
    return ;
}
void rebuild(){
    sorted_vertex.clear();
    sort_vertex();
    iterate_vertex();
    initial_request_for_ans();
}
void delete_edge(int u,int v){
    mp[make_pair(u,v)]=0;
    mp[make_pair(v,u)]=0;
    deg[u]--,deg[v]--;
    auto it=find(neighbor[u].begin(),neighbor[u].end(),v);
    if(it!=neighbor[u].end()){
        neighbor[u].erase(it);
    }
    it=find(neighbor[v].begin(),neighbor[v].end(),u);
    if(it!=neighbor[v].end()){
        neighbor[v].erase(it);
    }
    deg_gamma[u]=0,deg_gamma[v]=0;
    for(int w:neighbor[u]){
        deg_gamma[w]+=((deg[u]+1)>=(deg[w]+1)*ts)-((deg[u]+2)>=(deg[w]+1)*ts);
        deg_gamma[u]+=((deg[w]+1)>=(deg[u]+1)*ts);
    }
    for(int w:neighbor[v]){
        deg_gamma[w]+=((deg[v]+1)>=(deg[w]+1)*ts)-((deg[v]+2)>=(deg[w]+1)*ts);
        deg_gamma[v]+=((deg[w]+1)>=(deg[v]+1)*ts);
    }
    if(!use_bottom_k){
        for(int i=0;i<k;i++){
            in_case_empty(u,i,int(((long long)vec_a[i]*(long long)(v)+(long long)vec_b[i])%p),v);
            in_case_empty(v,i,int(((long long)vec_a[i]*(long long)(u)+(long long)vec_b[i])%p),u);
        }
    }
    else{
        remove_from_set(u,hash_bottom[v],v);
        remove_from_set(v,hash_bottom[u],u);
    }

    int i=0;
    while(i<sol_cnt){
        int num=sol[i].delete_how_much(u,v);
        if(num==3){
            if(i!=sol_cnt-1){
                sol[i]=sol[sol_cnt-1];
                sol_cnt--;
                continue;
            }
            else sol_cnt--;
        }
        i++;
    }
    return ;
}
void request_for_ans(){
    int ans=-1;
    deg_gamma_min=n;
    for(int i=0;i<sol_cnt;i++){
        if(sol[i].cur_sol_quasi>=alpha){
            if(ans==-1||sol[i].cur_sol_size>sol[ans].cur_sol_size){
                ans=i;
            }
        }
        deg_gamma_min=min(deg_gamma_min,sol[i].cur_sol_size);
    }
    if(ans==-1){
        for(int i=0;i<sol_cnt;i++){
            if(ans==-1||sol[i].cur_sol_quasi>=sol[ans].cur_sol_quasi){
                ans=i;
            }
        }
    }
    printf("%d %d %lf\n",sol[ans].cur_sol_vertex,sol[ans].cur_sol_size,sol[ans].cur_sol_quasi);
    return ;
}

void get_res(){
    clock_t begin = clock();
    o=0,e=0,sol_cnt=0;
    sol.resize(B);
    sol_vis.resize(n,0);
    sort_vertex();
    iterate_vertex();
    
    initial_request_for_ans();
    clock_t end = clock();
    double run_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << setprecision(4) << run_secs << endl;
    cout << "********************************" << endl;
    begin=clock();
    printf("%d\n",qccnt);
    qccnt=0,o=0,e=0,cnt=0;
    
    q=min(q,QLimited);
    int rebuild_cnt=0;
    int curtime=10;
    for(int i=1;i<=q;i++){
        int op,u,v;
        fin>>op>>u>>v;
        if(op==0) add_edge(u,v);
        if(op==1) delete_edge(u,v);
        if(sol_cnt==0||i%Batch==0){
            rebuild(),rebuild_cnt++;
        }
        request_for_ans();
    }
    printf("%d %d %d %d\n",o,e,rebuild_cnt,qccnt);
    end = clock();
    run_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << setprecision(4) << run_secs << endl;
    cout << "********************************" << endl;
    return ;
}

#endif
