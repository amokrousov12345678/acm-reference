// ************* Bridges, articulations points and blocks **************
vector<int> g[maxn];//IN g - undirected graph
int was[maxn];//MUST BE zeroed before each call on same component
int cur;
int num[maxn], low[maxn];

//find cutpoints in CONNECTED COMPONENT
int isCutpoint[maxn];//OUT if v is cutpoint
void dfsCutpoints(int v, int par = -1) {
    isCutpoint[v] = 0;
    was[v] = 1;
    num[v] = low[v] = cur++;
    int d = 0, maxlo = -1;
    for (auto& to: g[v]) {
        if (to == par) continue;
        if (was[to]) {
            low[v] = min(low[v], num[to]);
            continue;
        }
        dfsCutpoints(to, v); ++d;
        maxlo = max(maxlo, low[to]);
        low[v] = min(low[v], low[to]);
    }
    if (par >= 0) isCutpoint[v] = (maxlo >= num[v]);
    else isCutpoint[v] = (d > 1);
}
//find bridges in CONNECTED COMPONENT
//DOESN'T WORK WITH MULTIEDGES
vector<int> isBridge[maxn];//OUT isBridge[u][id] if edge g[u][id] is bridge
void dfsBridges(int v, int par = -1) {
	was[v] = 1;
	num[v] = low[v] = cur++;
	isBridge[v].clear(); isBridge[v].resize(g[v].size());
	for (int i = 0; i < g[v].size(); ++i) {
		isBridge[v][i] = 0;
		if (g[v][i] == par) continue;
		if (was[g[v][i]]) {
			low[v] = min(low[v], num[g[v][i]]);
			continue;
		}
		dfsBridges(g[v][i], v);
		low[v] = min(low[v], low[g[v][i]]);
		if (low[g[v][i]] > num[v]) isBridge[v][i] = 1;
	}
}
//TBD: write biconnected components in more simple way
int low[maxn], was[maxn], num[maxn], comp_num[maxn], cur;
int n;
vector<int> is_bridge[maxn], back[maxn], comp[maxn];
pair<int, int> stck[maxe];
int ssize, comp_cnt, is_art[maxn];

void add_biedge(int u, int v) {
	back[u].push_back(g[v].size());
	back[v].push_back(g[u].size());
	g[u].push_back(v);
	g[v].push_back(u);
	comp[u].push_back(-1);
	comp[v].push_back(-1);
	is_bridge[u].push_back(0);
	is_bridge[v].push_back(0);
}
void clear() {
	for (int i = 0; i < n; ++i) {
		g[i].clear();
		is_bridge[i].clear();
		back[i].clear();
		comp[i].clear();
	}
	memset(num, -1, sizeof(num));
	memset(was, 0, sizeof(was));
	cur = comp_cnt = ssize = 0;
}

// ***** Searching for biconnected components *****
void release(int u, int v) {
	while (ssize --> 0) {
		int i = stck[ssize].first, j = stck[ssize].second;
		comp[i][j] = comp[ g[i][j] ][ back[i][j] ] = comp_cnt;
		if ((i == u && g[i][j] == v) || (i == v && g[i][j] == u)) break;
	}
	comp_cnt++;
}
void dfs(int v, int par) {
	was[v] = 1;
	num[v] = low[v] = cur++;
	for (int i = 0; i < g[v].size(); ++i) {
		if (g[v][i] == par) continue;
		if (num[ g[v][i] ] < num[v]) stck[ssize++] = make_pair(v, i);
		if (was[ g[v][i] ]) {
			low[v] = min(low[v], num[ g[v][i] ]);
			continue;
		}
		dfs(g[v][i], v);
		low[v] = min(low[v], low[ g[v][i] ]);
		if (low[ g[v][i] ] >= num[v]) release(v, g[v][i]);
	}
}
