#include "NSF.h"


ifstream ifile;
void read_graph(string filename) {
    clock_t begin = clock();
    ifile.open(filename);
    if(ifile.fail())
    {
        cout << "### Error Open, File Name: " << filename << endl;
        return;
    }
    int u,v;
    ifile >> n >> m >> q;
    initial_variable();
    for (int i = 0; i < m; ++i) {
        ifile >> u >> v;
        degree[u]++;
        degree[v]++;
        neighbor[u].emplace_back(v);
        neighbor[v].emplace_back(u);
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

}

void InsertQuasiSizeAndDeleteOld(int node, int size){
    auto it = QuaSizeSet.find(make_pair(QuaSiz[node], node));
    if (it != QuaSizeSet.end()) {
        QuaSizeSet.erase(it);
    }
    QuaSiz[node] = size;
    QuaSizeSet.insert(make_pair(QuaSiz[node], node));
    return ;
}

void initial_variable() {
    QuaSiz.resize(n);
    QuaSizeSet.clear();
    degree.resize(n);
    neighbor.resize(n);
    kcore_size.resize(n);
    vis.resize(n);
    if (use_signature) {
        vec_a.resize(k);
        vec_b.resize(k);
        signature.resize(n);
    }
    for (int i = 0; i < k; i++) {
        vec_a[i] = rand() % (p - 1) + 1;
        vec_b[i] = rand();
    }
    for(int i=0;i<n;i++){
        QuaSiz[i]=0;
        vis[i]=0;
        degree[i]=0;
        QuaSizeSet.insert(make_pair(QuaSiz[i], i));
    }
}

void get_k_signatures(int node) {
    signature[node].clear();
    for (int index = 0; index < k; ++index) {
        int minhash = p;
        for (int j : neighbor[node]) {
            minhash = min(minhash, int(((long long)(vec_a[index]) * (long long)(j) + (long long)(vec_b[index])) % p));
        }
        minhash = min(minhash, int(((long long)(vec_a[index]) * (long long)(node) + (long long)(vec_b[index])) % p));
        signature[node].emplace_back(minhash);
    }
}

void sort_vertex() {
    for (int i = 0; i < n; ++i) {
        int num = core_num(i);
        sorted_vertex.emplace_back(pair<int,int>(num, i));
    }
    sort(sorted_vertex.begin(), sorted_vertex.end(), greater<pair<int,int>>());
}

void iterate_vertex() {
    for (int i = 0; i < sorted_vertex.size(); ++i) {
        auto p = sorted_vertex[i];
        int num = sift_num(p.second);
        InsertQuasiSizeAndDeleteOld(p.second, num);
        if (p.first > cur_sol_size || cur_sol_vertex == -1) {
            if (num > cur_sol_size || cur_sol_vertex == -1) {
                cur_sol_size = num;
                cur_sol_vertex = p.second;
            }
        }
    }
}

int core_num(int node) {
    int num = 0;
    for (int v : neighbor[node]) {
        if ((degree[v]+1) >= (degree[node]+1) * ts) num++;
    }
    return num;
}

double ct_score(int v1, int v2) {
    int res = 0;
    unordered_set<int> st;
    st.insert(v1);
    for (int num : neighbor[v1]) st.insert(num);
    for (int num : neighbor[v2]) if (st.count(num)) res++;
    if (st.count(v2)) res++;
    return (double)res / (double)(degree[v1] + 1);
}

double ct_score_minhash(int v1, int v2) {
    int res = 0;
    for (int i = 0; i < k; ++i) {
        if (signature[v1][i] == signature[v2][i]) res++;
    }
    double j_sim = (double)res / (double)k;
    return (double)(degree[v1] + degree[v2] + 2) * j_sim / (j_sim + 1) / (double)(degree[v1] + 1);
}

int sift_num(int node) {
    int res = 0;
    if (use_signature){
       if (vis[node] == 0) {
           get_k_signatures(node);
           vis[node] = 1;
       }
       for (int v : neighbor[node]) {
           if (vis[v] == 0) {
               get_k_signatures(v);
               vis[v] = 1;
           }
           if (ct_score_minhash(node, v) > ts) res++;
       }
    }
    else {
        for (int v : neighbor[node]) {
          if (ct_score(node, v) >= ts){
                res++;
            }
        }
    }
    double tmp = (double)res / (double)(degree[node]+1);
    if (tmp < b) return 0;
    return res + 1;
}
void add_edge(int u,int v){
    degree[u]++,degree[v]++;
    neighbor[u].emplace_back(v);
    neighbor[v].emplace_back(u);
    vis[u]=vis[v]=0;
    int num=local_detect(u);
    InsertQuasiSizeAndDeleteOld(u,num);
    num=local_detect(v);
    InsertQuasiSizeAndDeleteOld(v,num);
    return ;
}
void delete_edge(int u,int v){
    degree[u]--,degree[v]--;
    auto it=find(neighbor[u].begin(),neighbor[u].end(),v);
    if(it!=neighbor[u].end()){
        neighbor[u].erase(it);
    }
    it=find(neighbor[v].begin(),neighbor[v].end(),u);
    if(it!=neighbor[v].end()){
        neighbor[v].erase(it);
    }
    vis[u]=vis[v]=0;
    int num=local_detect(u);
    InsertQuasiSizeAndDeleteOld(u,num);
    num=local_detect(v);
    InsertQuasiSizeAndDeleteOld(v,num);
    return ;
}

void get_res(int flag = 0) {
    if(!flag){
        sorted_vertex.clear();
        for(int i=0;i<n;i++) vis[i]=0;
        cur_sol_size=0;
        cur_sol_edges=0;
        cur_sol.clear();
        cur_sol_vertex=-1;
        QuaSizeSet.clear();
        for(int i=0;i<n;i++){
            QuaSiz[i]=0;
            QuaSizeSet.insert(make_pair(QuaSiz[i], i));
        }
        sort_vertex();
        iterate_vertex();
    }
    cur_sol_size=0;
    cur_sol_edges=0;
    cur_sol.clear();
    cur_sol.emplace_back(cur_sol_vertex);
    for (int v : neighbor[cur_sol_vertex]) {
        if (use_signature) {
            if (ct_score_minhash(cur_sol_vertex, v) > ts) cur_sol.emplace_back(v);
        }
        else {
            if (ct_score(cur_sol_vertex, v) > ts) cur_sol.emplace_back(v);
        }
    }
    for (int i : cur_sol) {
        for (int i2 : neighbor[i]) {
            if(find(cur_sol.begin(), cur_sol.end(), i2) != cur_sol.end()) cur_sol_edges++;
        }
    }
    cur_sol_size = (int)cur_sol.size();
    if(cur_sol_size<=1) cur_sol_quasi=1.0;
    else cur_sol_quasi = double(cur_sol_edges) / double(cur_sol_size * (cur_sol_size-1));
}

int main(int argc, const char * argv[]) {
    if (argc < 3 || argc > 5){
        cout << "wrong input format!!" << endl;
        return 0;
    }
    if (argc == 5) {
        k = atoi(argv[4]);
        use_signature = true;
    }
    string filename;
    filename = argv[1];
    ts = atof(argv[2]);
    b = atof(argv[3]);
    Batch = 5000;

    read_graph(filename);

    clock_t begin = clock();
    get_res();
    clock_t end = clock();
    double run_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << setprecision(4) << run_secs << endl;
    cout << "********************************" << endl;
    begin=clock();
    
    printf("%d %d %lf\n",cur_sol_vertex,cur_sol_size,cur_sol_quasi);
    cout << cur_sol_size << " " << cur_sol_quasi << endl;
    q=q>100000?100000:q;
    
    int curtime=10;
    for(int i=1;i<=q;i++){
        int op,x,y;
        ifile>>op>>x>>y;
        if(i%Batch==0){
            if(op==0){
                add_edge(x,y);
            }
            else{
                delete_edge(x,y);
            }
            get_res();
        }
        else{
            if(op==0){
                add_edge(x,y);
                cur_sol_vertex=-1;
                int update_count=1;
                while((*QuaSizeSet.rbegin()).second!=cur_sol_vertex){
                    cur_sol_vertex=(*QuaSizeSet.rbegin()).second;
                    int num=local_detect(cur_sol_vertex);
                    InsertQuasiSizeAndDeleteOld(cur_sol_vertex,num);
                    update_count++;
                    if(update_count>20) break;
                }
                cur_sol_vertex=(*QuaSizeSet.rbegin()).second;
                get_res(1);
            }
            if(op==1){
                delete_edge(x,y);
                cur_sol_vertex=-1;
                int update_count=1;
                while((*QuaSizeSet.rbegin()).second!=cur_sol_vertex){
                    cur_sol_vertex=(*QuaSizeSet.rbegin()).second;
                    int num=local_detect(cur_sol_vertex);
                    InsertQuasiSizeAndDeleteOld(cur_sol_vertex,num);
                    update_count++;
                    if(update_count>20) break;
                }
                cur_sol_vertex=(*QuaSizeSet.rbegin()).second;
                get_res(1);
            }
        }
        cout << cur_sol_vertex << " " << cur_sol_size << " " << cur_sol_quasi << endl;
        if(i==curtime){
            end = clock();
            run_secs = double(end - begin) / CLOCKS_PER_SEC;
            cout << "Timing "  << setprecision(4) << run_secs << endl;
            curtime*=10;
        }
    }
    end = clock();
    run_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << setprecision(4) << run_secs << endl;
    cout << "********************************" << endl;
    return 0;
}