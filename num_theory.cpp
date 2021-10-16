// *************************** Number theory ***************************
int phi(int n) { //(NOT CHECKED)
    int res = n;
    for (int i = 2; i * i <= n; ++i) if (!(n % i)) {
        res /= n;
        res *= n - 1;
        for ( ; !(n % i); n /= i) ;
    }
    if (n > 1) { res /= n; res *= n - 1; }
    return res;
}
//Eratosphene sieve, O(nloglogn) time
char sieve[maxn]; //0 if prime, 1 - if not
for (int i = 2; i * i < maxn; i++) {
    if (sieve[i]) continue;
    for (int j = 2 * i; j < maxn; j += i) sieve[j] = 1;
}
//Linear sieve + factorizer, O(n) time
int lp[maxn];//last prime divisor (lp[n] = n => n is prime)
vector<ll> primes;
for (int i = 2; i < maxn; i++) {
    if (lp[i] == 0) {lp[i] = i; primes.push_back(i);}
    for (int j = 0; j < Sz(primes) && primes[j] <= lp[i] && i * primes[j] < maxn; ++j)
        lp[i * primes[j]] = primes[j];
}
//Calc mu(x) = (-1)^|primeDivisors| or 0 if contains prime^2. NEEDS lp FILLED BY SIEVE
int mu[maxn];
mu[1] = 1;
for (int i=2;i<maxn;i++) {
    ll pred = i/lp[i]; mu[i] = pred % lp[i] ? -mu[pred] : 0;
}
// discrete log PRIME MODULO (x, s.t. a^x=b(mod m)
int getLog(int a, int b, int m) { 
    int n = ((int) sqrt (m)) + 1;
    int an = 1; for (int i=0; i<n; ++i) an = mul(an, a, m);
    map<int,int> vals;
    for (int i=1, cur=an; i<=n; ++i) {
        if (!vals.count(cur)) vals[cur] = i;
        cur = mul(cur, an, m);
    }
    for (int i=0, cur=b; i<=n; ++i) {
        if (vals.count(cur)) {
            int ans = vals[cur] * n - i;
            if (ans < m) return ans;
        }
        cur = mul(cur,a,m);
    }
    assert(0); return -1;
}
int factmod(int n, int p) { //NOT CHECKED
    int res = 1; // factorial by modulo (with discarded pows of p)
    while (n > 1) {
        res = (res * ((n / p) % 2 ? p - 1 : 1)) % p;
        for (int i = 2; i <= n % p; ++i) res = (res * i) % p;
        n /= p;
    }
    return res % p;
}
//get x,y, s.t. a*x+b*y = gcd(a,b) (ret. value)
ll ext_gcd(ll a, ll b, ll& x, ll& y) {
    if (a == 0) { x = 0; y = 1; return b;}
    ll x1, y1;
    ll d = ext_gcd(b%a, a, x1, y1);
    x = y1 - (b / a) * x1; y = x1;
    return d;
}
//false - no solution, true - series X=x+k*dx, Y=y+k*dy, k in Z
bool solveDiofant(ll a, ll b, ll c, ll& x, ll& y, ll& dx, ll& dy) {
    ll _x, _y; ll gcd = ext_gcd(a,b,_x,_y);
    if (c % gcd) return false;
    x = _x*(c/gcd); y = _y*(c/gcd);
    dx = b/abs(gcd); dy = -a/abs(gcd);
    return true;
}