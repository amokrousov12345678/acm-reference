// ********************* Match in bipartite graph **********************
vector<int> g[maxn];//MUST BE undirected if edge (a,b) then should be (b,a)
int used[maxn];
int match[maxn];//match vert id (correct for verts from both parts)
int cc; //num of found edges

bool dfs(int v) {
    if (used[v] == cc) return false;
    used[v] = cc;
    for (auto& it: g[v]) {
        if (match[it]==-1 || dfs(match[it])) {
            match[v] = it; match[it] = v; return true;
        }
    }
    return false;
}

int kuhn(int n) {//try to match [1..n] (enough only one part, but may try ALL)
    memset(used, -1, n*sizeof(int));
    memset(match, -1, n*sizeof(int));
    cc = 0;
    for (int i=0;i<n;i++) {//may choose in other order (for example based on weight)
        if (match[i]!=-1) continue;
        if (dfs(i)) cc++;
    }
    return cc;
}

//TODO: add Hopcroft-Carp