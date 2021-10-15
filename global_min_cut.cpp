// *************************** Global mincut ***************************
// damages graph!, O(N^3)
int n;//IN n -verts cnt
cap_t g[maxn][maxn];//IN g[u][v] - weight (capacity) of (u,v). 
//CAUTION: g[u][v] = g[v][u]
vector<int> best_cut;//OUT verts of one part (remaining are in other part)
cap_t best_cost = INF;//OUT min sum weight of cutted edges between parts
//auxillary arrays:
vector<int> v[maxn];
int w[maxn];
bool exist[maxn], in_a[maxn];
void minCut() {
	best_cost = INFi; best_cut.clear(); 
	for (int i = 0; i < n; ++i) v[i].assign(1, i);
	memset(exist, true, sizeof(exist[0])*n);
	for (int ph = 0; ph < n - 1; ++ph) {
		memset(in_a, false, sizeof(in_a[0])*n);
		memset(w, 0, sizeof(w[0])*n);
		for (int it = 0, prev; it < n - ph; ++it) {
			int sel = -1;
			for (int i = 0; i < n; ++i)
				if (exist[i] && !in_a[i] && (sel == -1 || w[i] > w[sel])) sel = i;
			if (it == n - ph - 1) {
				if (w[sel] < best_cost) best_cost = w[sel],  best_cut = v[sel];
				v[prev].insert(v[prev].end(), v[sel].begin(), v[sel].end());
				for (int i = 0; i < n; ++i) g[prev][i] = g[i][prev] += g[sel][i]; 
				exist[sel] = false;
			}
			else {
				in_a[sel] = true;
				for (int i = 0; i < n; ++i) w[i] += g[sel][i];	
				prev = sel;
			}
		}
	}
}
