// Strongly-connected components
vector<int> g[maxn]; //IN g
int n;//IN n - vert count
int comp[maxn]; //OUT comp[x] - id of comp with vert x
//auxillary arrays
bool visit[maxn];
vector<int> revg[maxn];
vector<int> que, component;
int compCnt;

void dfs_fwd(int v) {
    visit[v] = true;
    for (int u : g[v]) if (!visit[u]) dfs_fwd(u);
    que.push_back(v);
}

void dfs_rev(int v) {
    visit[v] = true;
    for (int u : revg[v]) if (!visit[u]) dfs_rev(u);
    component.push_back(v);
}

void tarjan() {
	compCnt = 0;
    // add reverse edges
    for (int v = 0; v < n; ++v) revg[v].clear();
    for (int v = 0; v < n; ++v)
        for (int u : g[v]) 
            revg[u].push_back(v);
    // dfs forward
    memset(visit, 0, sizeof(visit[0])*n);
    que.clear();
    for (int v = 0; v < n; ++v) 
        if (!visit[v]) 
            dfs_fwd(v);
    reverse(que.begin(), que.end());
    // dfs backward
    memset(visit, 0, sizeof(visit[0])*n);
    for (int v : que) if (!visit[v]) {
        component.clear();
        dfs_rev(v);
        for (int u : component) comp[u] = compCnt;
        ++compCnt;
    }
}