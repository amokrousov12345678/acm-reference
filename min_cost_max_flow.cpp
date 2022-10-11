// ************************* Min-cost max-flow *************************
struct edge {
	int dest, back;//dest - vert to, back - id of reverse edge (actually id^1), 
	flow_t f, c; dist_t w;//OUT f - edge flow, IN c - capacity, w - edge cost
} edges[maxe];
int cnt, n;//IN cnt - cnt edges. Set to 0 before use on other graph, IN n - cnt verts
vector<int> g[maxn];//IN: g[v] - IDS of edges go from v

void add_edge(int u, int v, int c, dist_t w) {
	edges[cnt].dest = v; edges[cnt].back = cnt + 1;
	edges[cnt].f = 0; edges[cnt].c = c; edges[cnt].w = w; g[u].push_back(cnt);
	edges[cnt + 1].dest = u; edges[cnt + 1].back = cnt; edges[cnt + 1].f = 0;
	edges[cnt + 1].c = 0; edges[cnt + 1].w = -w; g[v].push_back(cnt + 1); cnt += 2;
}
int was[maxn], prevE[maxn];//was[v] - if v visited, prevE[v] - id of backedge to v
flow_t r[maxn]; dist_t dist[maxn], phi[maxn];//r[v] - min cap of path to v
//dist - sum cost on path, phi - johnson potential
bool dijk(int s, int t) {
	int i, j;
	for (i = 0; i < n; ++i) dist[i] = inf;
	memset(was, 0, sizeof(int) * n); memset(r, 0, sizeof(int) * n);
	dist[s] = 0; r[s] = inf;
	for (i = 0; i < n; ++i) {
		int mv = -1;
		for (j = 0; j < n; ++j) if (!was[j] && (mv == -1 || dist[mv] > dist[j])) mv = j;
		if (mv == -1 || dist[mv] == inf) break;
		was[mv] = 1;
		for (j = 0; j < g[mv].size(); ++j) {
			edge &e = edges[ g[mv][j] ];
			if (!was[e.dest] && e.f < e.c && dist[e.dest] > dist[mv] + e.w + phi[mv] - phi[e.dest]) {
				dist[e.dest] = dist[mv] + e.w + phi[mv] - phi[e.dest];
				r[e.dest] = min(r[mv], e.c - e.f); prevE[e.dest] = e.back;
			}
		}
	}
	return r[t] > 0;/*r[t] >= eps*/
}
dist_t aug(int s, int t) {
	flow_t rr = r[t]; dist_t ans = 0;
	while (s != t) {
		edge &e = edges[prevE[t]]; ans -= e.w * rr; e.f -= rr; edges[e.back].f += rr; t = e.dest;
	}
	return ans;
}

flow_t flow; dist_t cost;//OUT flow - max flow, cost - sum cost
void min_cost_max_flow(int s, int t) {
	cost = flow = 0; int i; for (i = 0; i < n; ++i) phi[i] = 0; 
	for (i = 0; i < cnt; ++i) edges[i].f = 0; 
	while (dijk(s, t)) { flow += r[t]; cost += aug(s, t); for (i = 0; i < n; ++i) phi[i] += dist[i];}
}