// ****** Hungarian algorithm (min cost assignment on n*m grid) ************************
int n, m; //IN n - row count, m - col count (MUST n<=m)
val_t c[maxn][maxn]; //IN: input matrix (1-based!!!!!)
int p[maxn]; // answer: p[i] for 1<=i<=m - id row assigned to col i (1-based!!!!)
//auxillary arrays:
int prv[maxn], nw[maxn];
val_t u[maxn], v[maxn], minv[maxn];//u,v, s.t. u[i]+v[j]<=c[i][j] and sum u + sum v -> max
val_t minCostAssignment() {//returns assignment cost
	memset(u, 0, sizeof(val_t) * (n + 1));
	memset(v, 0, sizeof(val_t) * (m + 1));
	memset(p, 0, sizeof(int) * (m + 1));
	memset(prv, 0, sizeof(int) * (m + 1));
	for (int i = 1; i <= n; ++i) {
		int j0 = 0, j1 = -1;
		p[0] = i;
		for (int j = 0; j <= m; ++j) minv[j] = INF;
		minv[0] = 0;
		memset(nw, 0, sizeof(int) * (m + 1));
		do {
			nw[j0] = 1;
			int i0 = p[j0];
			val_t del = INF;
			for (int j = 1; j <= m; ++j)
				if (nw[j] == 0) {
					val_t cc = c[i0][j] - u[i0] - v[j];
					if (cc < minv[j]) {
						minv[j] = cc; prv[j] = j0;
					}
					if (minv[j] < del) {
						del = minv[j]; j1 = j;
					}
				}
			for (int j = 0; j <= m; ++j)
				if (nw[j]) {
					u[p[j]] += del; v[j] -= del;
				} else
					minv[j] -= del;
			j0 = j1;
		} while (p[j0] != 0);
		do {
			j1 = prv[j0]; p[j0] = p[j1]; j0 = j1;
		} while (j0);
	}
	return -v[0];
}