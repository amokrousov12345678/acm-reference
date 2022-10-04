// ***************** Fast IO (noticeably faster than cin/cout)**********
//Printf/scanf almost no profit agains cin/cout with sync_with_stdio
#include <cstdio>
//Read
static const int buf_size = 4096;
inline int getChar() {
    static char buf[buf_size]; static int len = 0, pos = 0;
    if (pos == len) {pos = 0; len = fread(buf, 1, buf_size, stdin);}
    if (pos == len) return -1;
    return buf[pos++];
}
inline int readChar() {
    int c = getChar(); while (c <= 32) {c = getChar();} return c;
}
template <class T = int>
inline T readInt() {
    int s = 1, c = readChar(); T x = 0;
    if (c == '-') {s = -1; c = getChar();}
    while ('0' <= c && c <= '9') {x = x * 10 + c - '0'; c = getChar();}
    return s == 1 ? x : -x;
}
//Write
static int write_pos = 0; static char write_buf[buf_size];
inline void writeChar(int x) {
    if (write_pos == buf_size) {fwrite(write_buf, 1, buf_size, stdout); write_pos = 0;}
    write_buf[write_pos++] = x;
}
template <class T = int>
inline void writeInt(T x, char end = '\0') {
    if (x < 0) {writeChar('-'); x = -x;}
    char s[24]; int n = 0;
    while (x || !n) {s[n++] = '0' + x % 10; x /= 10;}
    while (n--) {writeChar(s[n]);} if (end) {writeChar(end);}
}
inline void writeWord(const char *s) { while (*s) writeChar(*s++);}
struct Flusher {
    ~Flusher() {if (write_pos) {fwrite(write_buf, 1, write_pos, stdout); write_pos = 0;}}
} flusher;//this MUST be in source to flush on app exit
// ***************** Setting the line up in the bitset *****************
inline const bool test(dword *arr, const int &x) {//getBit
	return (arr[x >> 5] >> (x & 31)) & 1;
}
inline void setBit(dword *arr, const int &x) {
	arr[x >> 5] |= (1 << (x & 31));
}
inline const dword getMask(const int &l, const int &r) {
	return (r < 32 ? (dword(1) << r) : dword(0)) - (dword(1) << l);
}
inline void setLine(dword *arr, int l, int r) {
	if (l > r) return;
	int li = (l >> 5);
	int ri = (r >> 5);
	if (li == ri) {
		arr[li] |= getMask(l - (li << 5), r + 1 - (ri << 5));
		return;
	}
	arr[li] |= getMask(l - (li << 5), 32);
	arr[ri] |= getMask(0, r + 1 - (ri << 5));
	for (int i = li + 1; i <= ri - 1; ++i) arr[i] = dword(-1);
}

// **************************** all submasks ***************************
for (int s = m; ;s = (s - 1) & m) {
    do_smth(s);
    if (!s) break;
}

//***************************** Geometry *******************
template<typename T>
struct PointT { //integer point
    T x, y;
    inline PointT() : x(0), y(0) {}
    inline PointT(T _x, T _y) : x(_x), y(_y) {}
    inline T len2() const { return x * x + y * y; }
    inline double len() const { return sqrt(len2()); }
    inline const PointT operator + (const PointT &b) const {
        return PointT(x + b.x, y + b.y);
    }
    inline const PointT operator - (const PointT &b) const {
        return PointT(x - b.x, y - b.y);
    }
    inline const PointT operator * (T b) const {
        return PointT(x * b, y * b);
    }
    inline bool operator == (const PointT &b) const {
        return x == b.x && y == b.y;
    }
    inline bool half() const { //true if in lower halfspace
        return (y < 0 || (y == 0 && x < 0));
    }
    inline bool operator < (const PointT &b) const { 
        bool th = half() ? 1 : 0;
        bool bh = b.half() ? 1 : 0;
        if (th ^ bh) return th < bh;
        auto pv = vect(*this, b);
        return pv > 0;
    }//integer compare by polar angle (counted from OX CCW)

    inline friend auto vect(const PointT &a, const PointT &b) {
        return a.x * b.y - a.y * b.x;
    }
    inline friend auto scal(const PointT &a, const PointT &b) {
        return a.x * b.x + a.y * b.y;
    }   
};
//Next functions work with Point=PointT<ll> (but you may use with doubles with some epsilon)
// decides whether <p> is inside the oriented CCW angle (from <a> to <b>) (including bounds)
// both <a> and <b> are nonzero
bool inSector(const Point &a, const Point &b, const Point &p) {
	auto vab = vect(a, b); if (!vab && (scal(a, b) > 0)) if (scal(a, p) < 0) return false;
	if (vab >= 0) return (vect(a, p) >= 0 && vect(p, b) >= 0);
	else		  return (vect(a, p) >= 0 || vect(p, b) >= 0);
}
template<typename T> inline bool Cross1DSegs(T l1, T r1, T l2, T r2) {//1 dimensional
	if (l1 > r1) swap(l1, r1); if (l2 > r2) swap(l2, r2); return !(l1 > r2 || r1 < l2);
}
// (a1,b1 - one seg, a2,b2 - other seg)
bool CrossSegs(const Point &a1, const Point &b1, const Point &a2, const Point &b2) {
    auto a11 = b1.x - a1.x; auto a12 = a2.x - b2.x;
    auto a21 = b1.y - a1.y; auto a22 = a2.y - b2.y;
    auto xb1 = a2.x - a1.x; auto xb2 = a2.y - a1.y;
    auto det  = a11 * a22 - a12 * a21; auto detu = xb1 * a22 - a12 * xb2;
    auto detv = a11 * xb2 - xb1 * a21;
    if (det == 0) {
        if (detu || detv) return false;
        if (a1.x != b1.x || a2.x != b2.x) return Cross1DSegs(a1.x, b1.x, a2.x, b2.x);
        if (a1.y != b1.y || a2.y != b2.y) return Cross1DSegs(a1.y, b1.y, a2.y, b2.y);
        return (a1 == a2);
    }
    if (det < 0) {det = -det; detu = -detu; detv = -detv; }
    return (detu >= 0 && detu <= det && detv >= 0 && detv <= det);
}
// determines whether point <p> lies on the closed segment <a>-<b>
inline bool onSeg(const Point &p, const Point &a, const Point &b) {
	if (vect(p - a, b - a)) return false;   if (a == b) return p == a;
	return (scal(p - a, b - a) >= 0 && scal(p - b, a - b) >= 0);
}

// determines whether a point <p> is inside the triangle <a>-<b>-<c>
// does not work for triangles with zero area (NOT CHECKED)
bool inTriangle(const Point &p, const Point &a, const Point &b, const Point &c, bool strict = false) {
	auto tv = abs(vect(c - a, b - a)); auto t = abs(vect(p - a, b - a)); 
    if (strict && !t) return false;
	tv -= t; t = abs(vect(p - b, c - b)); if (strict && !t) return false;
	tv -= t; t = abs(vect(p - c, a - c)); if (strict && !t) return false;
	tv -= t; return tv >= 0;
}

// determines whether point <p> lies inside polygon <arr[0], ..., arr[n]>
// no self-crossings! But may be non-convex. May have equal points	arr[0] = arr[n]
bool inPolygon(const Point &p, int n, const Point *arr, bool strict = false) {
    //	if lies on the border
	for (int i = 0; i < n; ++i) if (onSeg(p, arr[i], arr[i + 1])) return !strict;
	Point spot = p + Point(15013, 15017); // BIG PRIMES: 1061109589 / 1061109601
	int cnt = 0;
	for (int i = 0; i < n; ++i) if (CrossSegs(spot, p, arr[i], arr[i + 1])) cnt ^= 1;
	return bool(cnt);
}
// Graham's convex hull
// changes order of points! no equal points allowed! must have nonzero area
Point hctr;
bool cmpHull(const Point &a, const Point &b) {
    Point ad = a - hctr; Point bd = b - hctr;
    auto tv = vect(ad, bd); if (tv) return tv > 0; //abs(tv) > eps for real
    return ad.len2() < bd.len2();
}
//out k - size of convex hull in res, points CCW in res
void ConvexHull(int n, Point *arr, int &k, Point *res, bool strict = true) {
	if (n<=2) {k = n; memcpy(res, arr, sizeof(Point)*n); return;};
    int i, best = 0;
    for (i = 1; i<n; i++) {
        if (arr[i].x < arr[best].x) best = i;
        if (arr[i].x == arr[best].x && arr[i].y < arr[best].y) best = i;
    }
    std::swap(arr[best], arr[0]); hctr = arr[0]; std::sort(arr+1, arr+n, cmpHull);
    k = 0; res[k++] = arr[0]; res[k++] = arr[1];
    for (i = 2; i<n; i++) {
        if (strict) while (k>=2 && vect(res[k-1]-res[k-2], arr[i]-res[k-2]) <= 0) k--;
        if (!strict) while (k>=2 && vect(res[k-1]-res[k-2], arr[i]-res[k-2]) < 0) k--;
        res[k++] = arr[i];
    }
    if (!strict) {
        k--;
        for (i = n-1; i>0; i--) {
                res[k++] = arr[i];
            if (vect(arr[i]-arr[0], arr[i-1]-arr[0]) != 0) break;
        }
    }
}

//NEXT FUNCTIONS work ONLY with Point=PointT<double>
//Line crosses Line (a1,b1 - one line, a2,b2 - other line) (infinite)
bool CrossLineLine(const Point &a1, const Point &b1, const Point &a2, const Point &b2, Point &res) {
	auto a11 = b1.x - a1.x;
	auto a12 = a2.x - b2.x;
	auto a21 = b1.y - a1.y;
	auto a22 = a2.y - b2.y;
	auto xb1 = a2.x - a1.x;
	auto xb2 = a2.y - a1.y;
	auto det  = a11 * a22 - a12 * a21;
	auto detu = xb1 * a22 - a12 * xb2;
	auto detv = a11 * xb2 - xb1 * a21;
	if (abs(det) < eps) return false; //parallel (detu!=0) or coincide (detu=0)	
	detu /= det;
	detv /= det;
	Point c = a1 + (b1 - a1) * detu;
	res = c;
	return true;
}
//Even in case of touch, later functions return 2 points, you should check if they coincide
//Line (la, lb) crosses Circle (center cc, radius cr) (infinite)
bool CrossLineCircle(const Point &la, const Point &lb, const Point &cc, real_t cr, Point &res1, Point &res2) {
	Point st = la - cc; Point dir = lb - la;
	auto qa = scal(dir, dir); auto qb = 2.0 * scal(st, dir); 
    auto qc = scal(st, st) - cr * cr; auto qd = qb * qb - 4.0 * qa * qc;
	if (qd < -eps) return false; 
	if (qd < 0.0) qd = 0.0; qd = sqrt(qd); 
    auto x1 = (-qb - qd) / (2.0 * qa); auto x2 = (-qb + qd) / (2.0 * qa);
	res1 = la + dir * x1; res2 = la + dir * x2; return true;
}   
//Circle (center c1, radius r1) crosses Circle (center c2, radius r2)
bool CrossCircleCircle(const Point &c1, real_t r1, const Point &c2, real_t r2, Point &res1, Point &res2) {
	auto la = 2.0 * (c2.x - c1.x);
	auto lb = 2.0 * (c2.y - c1.y);
	auto lc = sqr(c1.x) - sqr(c2.x) + sqr(c1.y) - sqr(c2.y) + sqr(r2) - sqr(r1);
	if (la * la + lb * lb < eps) return false; //circle or coincide or no intersection (based on radius)
	Point a, b;
	if (abs(la) > abs(lb)) { //Logic can be used to build 2 points based on normal line equation
		a = Point(-lc / la, 0.0); b = Point(-(lb + lc) / la, 1.0);
	}
	else { a = Point(0.0, -lc / lb); b = Point(1.0, -(lc + la) / lb);}
	return CrossLineCircle(a, b, c1, r1, res1, res2);
}
//*******************************************HASHES******************************
//For set: sort items and get polyhash of this sequence
//For rooted tree: for leaf - some const, for node - hash of sorted child hashes sequence
//CAUTION: polyhash don't work here, use [sum(i=0 to k)((a_i)^2+a_i*p^i+3)] mod 2^64
//Or store in hashset polyhash of vector<H> for children. Hash of some subtree - id in set (1-started)
//Sometimes useful to use non-poly hash is you don't need to extract substr
namespace HPar {//when change SZ MUST copypaste logic for other base(s) EVERYWHERE
    const int SZ = 1;//Primitives: (3(27)(243) by 998244353, 5(125) by 1e9+7, 13 by 1e9+9)
    int mods[SZ], base[SZ], invBase[SZ], basePows[SZ][maxn], invBasePows[SZ][maxn];
    void init() {
        mods[0] = 998244353; base[0] = 27; invBase[0] = inv(base[0], mods[0]);
        basePows[0][0] = 1; invBasePows[0][0] = 1;
        for (int i=1;i<maxn;i++) {
            basePows[0][i] = mul(basePows[0][i-1], base[0], mods[0]);
            invBasePows[0][i] = mul(invBasePows[0][i-1], invBase[0], mods[0]);
        }
    }
}
struct H {
    std::array<int, HPar::SZ> data{0};
    void append(int pos, char ch) {//pos - RELATIVE in hash, coef
        data[0] = add(data[0], mul(HPar::basePows[0][pos], numc(ch), HPar::mods[0]), HPar::mods[0]);
    }
    void erase0th(char ch) {//ch - CH on this pos, all other "chars" shift left
        data[0] = sub(data[0], numc(ch), HPar::mods[0]);
        data[0] = mul(data[0], HPar::invBase[0], HPar::mods[0]);
    }
    friend H operator-(const H& l, const H& r) {
        H res; res.data[0] = sub(l.data[0], r.data[0], HPar::mods[0]); return res;
    }
    void shiftLeft(int cnt) {//indexing [l;r] -> [l-cnt;r-cnt] (if [l-cnt;l) initially empty!!)
        data[0] = mul(data[0], HPar::invBasePows[0][cnt], HPar::mods[0]);
    }
    friend bool operator==(const H& l, const H& r) {return l.data==r.data;}
};
namespace std { template<> struct hash<H> {
        std::size_t operator()(H const& h) const noexcept {return h.data[0]/* ^ h.data[1]*/;}
};};
vector<H> getPrefHashes(const string& s) {//res[i] - hash of pref [0..i]
    vector<H> res; H acc;
    for (int i=0;i<Sz(s);i++) {acc.append(i, s[i]); res.push_back(acc);} return res;
}
H getSubstrHash(const vector<H>& prefHashes, int l, int r) {//get [l;r] hash
    assert(l<=r && r<Sz(prefHashes));
    H res = prefHashes[r]; if (l) res = res-prefHashes[l-1];
    res.shiftLeft(l); return res;
}//Sometimes binsearch in array faster than seach in unordered_set
//In add,sub,mul avoid mod ops as possible (in mul guard by if works faster in this case)
//Window processing may be faster then prefix precalc and extract range
//*****************************Euler tour*******************
vector<int> path;//MUST clear before use. Exists if number of odd verts 0 (or 2, then patch by edge)
void dfs(int v) {//MUST reverse printed cycle when orinted edges
    while (!g[v].empty()) {
        int to = *g[v].begin(); g[v].erase(to);/*g[to].erase(v); (for bidir edges)*/dfs(to);
    }
    path.push_back(v);
}
// ******************************* Gauss *********************
int n, m;//IN n - rows, m - cols
real_t matr[maxn][maxn];//IN/OUT matr[N][M]. For Ax=B equation put B in m+1 col
//if without B, replace <=m to <m (and u=m to u=m-1) (4 times total)
int r;//OUT: count of non-zero rows after diag (rank)
int adr[maxn];//id of non-null column in row
bool used[maxn];//OUT: if column non-empty (i.e. var not free)
//if calc det, remember where do swap or normalization row
//if integer modulo, remove EPS, divide modulo. If modulo=2, use bit ops
//to invert write in matr 2n*n (A|E), and diag. Result will be (E|A^(-1))
void Gauss() {//diagonalize matrix
	int i, j, u;
	r = 0; //first row of remaining part
	memset(used, 0, sizeof(used));
	for (i = 0; i <= m; ++i) {
		int best = -1;
		for (j = r; j < n; ++j) if (best < 0 || abs(matr[j][i]) > abs(matr[best][i])) best = j;
		if (best < 0) break;//no rows left, stop
		for (u = 0; u <= m; ++u) swap(matr[best][u], matr[r][u]);
		if (abs(matr[r][i]) < EPS) continue;//current column is zero, skip it
		for (u = m; u >= i; --u) matr[r][u] /= matr[r][i]; //norm row
		for (j = 0; j < n; ++j) if (j != r) {
				real_t coef = matr[j][i];
				for (u = i; u <= m; ++u) matr[j][u] -= coef * matr[r][u];
			}//nullify all poses in column except of diagonal
		//may ignore upper/lower to get trigonal matrix little bit faster
		used[i] = true;
		adr[r++] = i;
	}
}
real_t sol[maxn];//var values which is one of solutions (if exists any)
bool GetSolution() {//uses matr which should be upper trigonal
	int i, j;
	memset(sol, 0, sizeof(sol));
	if (used[m]) return false;//incompatible system
	sol[m] = -1.0;    //MUST BE SO!
	for (i = 0; i < m; ++i)
		if (!used[i]) sol[i] = rand() / 32768.0; //free variables (any vals)
	for (i = r - 1; i >= 0; --i)
		for (j = adr[i] + 1; j <= n; ++j)
			sol[adr[i]] -= sol[j] * matr[i][j];
	return true;
}
// ******************************** LCA ********************************
// LCA of two nodes: <a> and <b>, needs (N*logN)*sizeof(int) memory
int n;//IN graph sz
int h[maxn]; //IN distance from root
int par[LOGmaxn/*+1*/][maxn];//IN par[0][x] - par of vert (par[0][root]=root
void LCAInit() {//OUT par - lifts
	for (int i = 1; i < LOGmaxn; ++i) {
		for (int j = 0; j < n; ++j) {
			par[i][j] = par[i - 1][par[i - 1][j]];
		}
	}
}
int LCA(int a, int b) {
	if (h[a] < h[b]) swap(a, b);
	int delta = h[a]-h[b];
	for (int i = LOGmaxn - 1; i >= 0; --i)
		if ((delta >> i) & 1) a = par[i][a];
	for (int i = LOGmaxn - 1; i >= 0; --i) if (par[i][a] != par[i][b]) {
		a = par[i][a];
		b = par[i][b];
	}
	if (a != b) a = par[0][a];
	return a;
}
//Rq time O(log n). For fast write euler tour of tree  (write vertex on enter and 
//after each child) and do RMQ via sparse table (look to suffmass.cpp)
// ******************************* DSU **********************************
int parent[maxn], sz[maxn];//parent link and component sz
void make_set(int v) {//MUST be called for each v before usage
    parent[v] = v; sz[v] = 1;
}
int find_set(int v) {//find_set(a)=find_set(b) if a and b in one set
    if (v == parent[v]) return v;
    return parent[v] = find_set (parent[v]);
}
void union_sets(int a, int b) {//join sets, which contain a and b
    a = find_set (a); b = find_set (b);
    if (a != b) {
        if (sz[a] < sz[b]) swap (a, b);
        parent[b] = a;
        sz[a] += sz[b];
    }
}
//given arr[0..2^bits-1]. Calc res[i] = sum(j submask of i)arr[j]. Stored in same arr. 
//May be used ^=, &=, |= instead of +- (any commutative associative ops)
void sumOnSubsets(T* arr, int bits) {
    int sz = 1<<bits;
    for (int j=bits-1;j>=0;j--) {
        int shift = 1<<(j+1); int len = 1<<j;
        for (int i=0;i<sz;i+=shift) for (int k=i;k<i+len;k++) arr[k+len] += arr[k];
    }
}
//CONVEX HULL TRICK: calculate min/max (ki*x+bi) on set lines (ki,bi) with convex hull idea
//for min - MUST add in decreasing order by k, for max - increasing (draw to understand)
template<bool rToL = false> //to add in reverse order, pass rToL = true
struct CHT {
    vector<line> lines;
    vector<ll> ints; //ord numbers for lower bound
    void add(const line& l) {
        while (lines.size()>=2) {
            double xl = inters(lines[Sz(lines)-2], lines[Sz(lines)-1]);
            double xr = inters(lines[Sz(lines)-1], l);
            if ((xl > xr && !rToL) || (xl < xr && rToL)) lines.pop_back(); else break;
        }
        ints.resize(lines.size()); ints.push_back(Sz(ints)); lines.push_back(l);
    }
    ll rq(ll x) {
        assert(!lines.empty());
        auto it = *lower_bound(ints.begin(), ints.end()-1, x,[this](ll id, ll x) {
            if (!rToL) return inters(lines[id], lines[id+1]) < x;
            else return inters(lines[id], lines[id+1]) > x;
        });
        return lines[it].f(x);
    }
};
// **************************** maxflow:lift ****************************
//IN n - verts cnt, c[u][v] - (u,v) capactity, s - source, t - target
int n, s, t; 
flow_t c[maxn][maxn];//for val_t with float point, use logic with eps
//O(V^2*sqrt(E)), in practice often slower than Dinitz on arbitrary graph
flow_t f[maxn][maxn]; //OUT f[u][v] - (u,v) flow
flow_t e[maxn]; int h[maxn], maxh[maxn];
	
flow_t pushRelabel() {//returns flow amount
	for (int i=0;i<n;i++) for (int j=0;j<n;j++) f[i][j] = 0;
	for (int i = 0; i < n; ++i) {
		f[s][i] = c[s][i]; f[i][s] = -f[s][i];
		e[i] = c[s][i]; h[i] = 0; maxh[i] = 0;
	}
	h[s] = n - 1; 
	int sz = 0;
	for (;;) {
		if (!sz)
			for (int i = 0; i < n; ++i) if (i != s && i != t && e[i] > 0) {
				if (sz && h[i] > h[maxh[0]]) sz = 0;
				if (!sz || h[i] == h[maxh[0]]) maxh[sz++] = i;
			}
		if (!sz) break;
		while (sz) {
			int i = maxh[sz - 1];
			bool pushed = false;
			for (int j = 0; j < n && e[i]; ++j)
				if (c[i][j] - f[i][j] > 0 && h[i] == h[j] + 1) {//abs(c-f)> eps
					pushed = true;
					flow_t addf = min(c[i][j] - f[i][j], e[i]);
					f[i][j] += addf, f[j][i] -= addf;
					e[i] -= addf; e[j] += addf;
					if (e[i] == 0) --sz;
				}
			if (!pushed) {
				h[i] = 1e9;
				for (int j = 0; j < n; ++j) //abs(c-f)> eps
					if (c[i][j] - f[i][j] > 0 && h[j] + 1 < h[i]) h[i] = h[j] + 1;
				if (h[i] > h[maxh[0]]) {
					sz = 0;
					break;
				}
			}
		}
	}
	flow_t flow = 0; //result flow value
	for (int i=0;i<n;i++)/* if (c[i][t])*/ flow += f[i][t]; 
	return flow;
}
// *************************** maxflow:Dinits ***************************
//IN n - verts cnt, c[u][v] - (u,v) capactity, s - source, t - target
int n, s, t; 
flow_t c[maxn][maxn];//for val_t with float point, use logic with eps
//O(V^2*E), O(E*sqrtV) on unit network (matching)
flow_t minPushed = 1;//INTERNAL min flow we can push through edge now. MUST BE >0
//(1 for int, eps for real, or changed automatically if dinicScale called)
flow_t f[maxn][maxn]; //OUT (u,v) flow
int d[maxn], ptr[maxn], q[maxn];
bool bfs() {
	int qh=0, qt=0;
	q[qt++] = s;
	memset(d, -1, n * sizeof(d[0]));
	d[s] = 0;
	while (qh < qt) {
		int v = q[qh++];
		for (int to = 0; to < n; ++to) if (d[to] == -1 && f[v][to] + minPushed <= c[v][to]) {
			q[qt++] = to;
			d[to] = d[v] + 1;
		}
	}
	return d[t] != -1;
}
flow_t dfs (int v, flow_t flow) {
	if (flow < minPushed) return 0;
	if (v == t) return flow;
	for (int &to = ptr[v]; to < n; ++to) {
		if (d[to] != d[v] + 1)  continue;
		flow_t pushed = dfs(to, min(flow, c[v][to] - f[v][to]));
		if (pushed >= minPushed) {
			f[v][to] += pushed;
			f[to][v] -= pushed;
			return pushed;
		}
	}
	return 0;
}
flow_t dinic(bool doClean = true) {//returns flow amount. doClean - should we reset f[i][j]
	if (doClean) for (int i=0;i<n;i++) for (int j=0;j<n;j++) f[i][j] = 0;
	flow_t flow = 0;
	for (;;) {
		if (!bfs()) break;
		memset(ptr, 0, n * sizeof ptr[0]);
		while (flow_t pushed = dfs (s, INF)) flow += pushed;
	}
	return flow;
}
flow_t dinicScale() {//returns flow amount, uses SCALING IDEA, modifies minPushed
	minPushed = INF; flow_t flow = 0;
	while (minPushed) { /*minPushed > eps*/
		flow += dinic(false); minPushed /= 2;
	}
	return flow;
}