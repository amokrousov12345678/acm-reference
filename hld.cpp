int par[maxn]; //parent
int sz_sub[maxn]; //subtree sz
int mx[maxn]; //"heavy" child (actually non-light) or "-1"
int pos[maxn]; //pos of vertex in segtree
int root[maxn]; //top-most vert of heavy chain
int depth[maxn]; //depth from root
void dfs_all(int v, int d = 0) {//first this
	depth[v] = d;
	sz_sub[v] = 1;
	int ind = -1, mai = -1;
	for (auto it : g[v]) {
		if (par[v] != it ) {
			par[it] = v;
			dfs_all(it, d + 1);
			sz_sub[v] += sz_sub[it];
			if (sz_sub[it] > mai) {
				mai = sz_sub[it];
				ind = it;
			}
		}
	}
	mx[v] = ind;	
}
void build_array () {//then this
	int cpos = 0;
	for (int i = 0; i < maxn; ++i) { //sz_sub[i], cause we MUST ignore not-inited verts with zeros in arrays
		if (sz_sub[i] && (par[i] == -1 || mx[par[i]] != i)) { 
			for (int v = i; v != -1; v = mx[v]) {
				pos[v] = cpos++;
				root[v] = i;
			}
		}
	}
}
int get(int u, int v) {
	int ans = 0;
	while (root[u] != root[v]) {
		if (depth[root[u]] < depth[root[v]]) swap(u, v);
		ans = max(ans, rmq.query(pos[root[u]], pos[u])) ;
		u = par[root[u]];
	}
	if (depth[u] < depth[v]) swap(u, v) ;
	ans = op(ans, rmq.query(pos[v], pos[u]));
	return ans ;
}
