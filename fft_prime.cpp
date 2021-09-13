const int mod = 998244353;//7340033;
const int root = 15311432 /*3^(7*17)*/ //5;
const int root_1 = //4404020;
const int root_pw = //1<<20;
 
//n MUST BE 2^k
void fft (int* a, int n, bool invert) {
	int n = (int) a.size();
 
	for (int i=1, j=0; i<n; ++i) {
		int bit = n >> 1;
		for (; j>=bit; bit>>=1)
			j -= bit;
		j += bit;
		if (i < j)
			swap (a[i], a[j]);
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

//inv series: A(x)*Bn(x) = 1 + x^n*C(x): B1 = A(x)[0]^-1, B2k = Bk*(2-A*Bk) (only 2k members)
//division: deg A = m, degB = n. A = BQ + R, Br, Ar, Qr - polys with reverse coeff order,
//D - (m+1-n) first of Br. Then Qr = D*A (m-n+1 coeffs, need to reverse)