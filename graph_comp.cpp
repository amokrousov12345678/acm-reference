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
//Biconn. component of verts - set of verts, each pair have at least 2 EDGE independent paths
//Can be found via dfs, when traverse bridge - enter new component
//Biconn. component of edges - set of edges, each pair have at least 2 VERT independent paths between ends
//Found by algo below
vector<int> compId[maxn];//OUT compId[u][id] - id of biconn component of edge g[u][id]

vector<int> back[maxn];
pair<int, int> stck[maxe];
int ssize = 0;//MUST be set to 0 before call
void release(int u, int v) {
	while (ssize-- > 0) {
		int i = stck[ssize].first, j = stck[ssize].second;
		comp[i][j] = comp[g[i][j]][back[i][j]] = comp_cnt;
		if ((i == u && g[i][j] == v) || (i == v && g[i][j] == u)) break;
	}
	comp_cnt++;
}
void dfs(int v, int par) {//Biconn. component of edges in CONNECTED COMPONENT
	compId[v].clear(); compId[v].assign(g[v].size(), -1);
	was[v] = 1;
	num[v] = low[v] = cur++;
	for (int i = 0; i < g[v].size(); ++i) {
		if (g[v][i] == par) continue;
		if (num[g[v][i]] < num[v]) stck[ssize++] = make_pair(v, i);
		if (was[g[v][i]]) {
			low[v] = min(low[v], num[g[v][i]]);
			continue;
		}
		dfs(g[v][i], v);
		low[v] = min(low[v], low[g[v][i]]);
		if (low[g[v][i]] >= num[v]) release(v, g[v][i]);
	}
}