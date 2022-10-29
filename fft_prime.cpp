//FFT<7340033, 78125/*5^7*/, 1<<20>; FFT<469762049, 2187/*3^7*/, 1<<26>;
//FFT<897581057, 872686320/*3^107*/, 1<<23>; FFT<880803841, 98722167/*13^(3*5*7)*/, 1<<23>;
//FFT<1004535809, 702606812/*3^479*/, 1<<21>; FFT<1012924417, 673144645/*5^(3*7*23)*/, 1<<21>;
template<int mod=998244353, int root=15311432/*3^(7*17)*/, int root_pw=1<<23> struct FFT {
    const static int mod_ = mod; const int root_1 = inv(root, mod);
    void fft (int* a, int n, bool invert) {//n MUST BE 2^k
        assert(n<=root_pw);
        for (int i=1, j=0; i<n; ++i) {
            int bit = n >> 1;
            for (; j>=bit; bit>>=1)
                j -= bit;
            j += bit;
            if (i < j) swap (a[i], a[j]);
        }
        for (int len=2; len<=n; len<<=1) {
            int wlen = invert ? root_1 : root;
            for (int i=len; i<root_pw; i<<=1)
                wlen = int (wlen * 1ll * wlen % mod);
            for (int i=0; i<n; i+=len) {
                int w = 1;
                for (int j=0; j<len/2; ++j) {
                    int u = a[i+j],  v = int (a[i+j+len/2] * 1ll * w % mod);
                    a[i+j] = u+v < mod ? u+v : u+v-mod;
                    a[i+j+len/2] = u-v >= 0 ? u-v : u-v+mod;
                    w = int (w * 1ll * wlen % mod);
                }
            }
        }
        if (invert) {
            int nrev = inv(n, mod);
            for (int i=0; i<n; ++i)
                a[i] = int (a[i] * 1ll * nrev % mod);
        }
    }
    vector<int> mul(vector<int> a, vector<int> b) {
        int sumLen = Sz(a) + Sz(b) - 1;
        int log2 = 0; while ((1<<log2)<sumLen) log2++;
        a.resize(1<<log2, 0); b.resize(1<<log2, 0);
        fft(a.data(), Sz(a), false); fft(b.data(), Sz(b), false);
        for (int i=0;i<(1<<log2);i++) a[i] = ::mul(a[i], b[i], mod);
        fft(a.data(), Sz(a), true);
        return a;
    }
};
//FWHT: (a b) => (a+b a-b) (for inv should div by 2)	
//inv series: A(x)*Bn(x) = 1 + x^n*C(x): B1 = A(x)[0]^-1, B2k = Bk*(2-A*Bk) 
//from A taken 2k summands, mul in 4k window and select only 2k summands from result
//if multidimensial expand by first, than by second, etc
//(by completed dimension select only N summands is enough)
//division: deg A = m, degB = n. A = BQ + R, Br, Ar, Qr - polys with reverse coeff order,
//D - (m+1-n) first of Br. Then Qr = D*A (m-n+1 coeffs, need to reverse)