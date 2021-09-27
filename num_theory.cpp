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

int get_log(int a, int b, int m) { // discrete log (x, s.t. a^x=b(mod m)
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
    int res = 1; // factorial by modulo (with thrown pows of p)
    while (n > 1) {
        res = (res * ((n / p) % 2 ? p - 1 : 1)) % p;
        for (int i = 2; i <= n % p; ++i) res = (res * i) % p;
        n /= p;
    }
    return res % p;
}
