// finds matching, not maximal with N / MOD probability (2e-6 on n=500);
int mat[maxn][maxn], a[maxn][maxn], idx[maxn], iidx[maxn], has[maxn];//auxillary arrays
mt19937 rnd(chrono::steady_clock::now().time_since_epoch().count());
uniform_int_distribution<> distr(1, mod - 1);
//n - number of verts, ed - edges of graph (as pair of ends). Return edges in matching in same form
vector<pair<int, int>> matching(int n, vector<pair<int, int>>& ed)
{
    for (int i=0;i<n;i++) {idx[i]=i; iidx[i]=i; has[i]=0;}
    for (auto e : ed)
    {
        int u = idx[e.first]; int v = idx[e.second];
        if (u != -1 && v != -1) {int r = distr(rnd); mat[u][v] = r; mat[v][u] = mod - r;}
    }
    //inv considers matr as n*n. Return rank of matrix. If rank=n, set a=inv(matr).
    int r = inv(n); int m = 2 * n - r; assert(r % 2 == 0);// May be done by gaussing (matr|E)
    if (m != n)
    {
        int t = 0;
        do {
            memset(tmp, 0, sizeof(tmp));
            for (int i = 0; i < n; i++) for (int j = n; j < m; j++) {
                int r = distr(rnd);
                mat[i][j] = r, mat[j][i] = mod - r;
            }
        } while (inv(m) != m);
    }
    vector<pair<int, int>> ret;
    int fi, fj;
    for (int it = 0; it < m / 2; it++) {
        bool found = 0;
        for (int i = 0; !found && i < m; i++)
            if (!has[i]) for (int j = i + 1; j < m; j++)
                    if (a[i][j] && mat[i][j]) {fi = i; fj = j; found = 1; break;}
        assert(found); if (fj < n) ret.push_back({ iidx[fi], iidx[fj] });
        has[fi] = has[fj] = 1;
        for (int sw = 0; sw < 2; sw++)
        {
            ll v = binPow(a[fi][fj], mod - 2);//normal binPow modulo mod
            for (int i = 0; i < m; i++) if (!has[i] && a[i][fj])
            {
                ll b = a[i][fj] * v % mod;
                for (int j = 0; j < m; j++){
                    a[i][j] -= a[fi][j] * b % mod; if (a[i][j] < 0) a[i][j] += mod;
                }
            }
            swap(fi, fj);
        }
    }
    return ret;
}