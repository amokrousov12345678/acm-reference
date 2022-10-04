// ************************* String functions **************************
inline int numc(char c) { return c <= 'Z' ? c - 'A' : c - 'a' + 26; }
void prefix_function(const char *s, int *pi) {
    //p[i] = max len of non-trivial (len!=i+1) prefix = suffix of S[0..i] (or 0 if such not exist)
	int n = strlen(s);
	for (int i = 1; i < n; ++i) {
		int j = pi[i - 1];
		while (j > 0 && s[i] != s[j]) j = pi[j - 1];
		if (s[i] == s[j]) ++j;
		pi[i] = j;
	}
}
void z_function(const char *st, int *z) {//Z[i] = max common prefix of S[0..n-1] and S[i..n-1]
	int i, j = 0, r = 0, l;
	z[0] = l = strlen(st);
    //z[l] = '#';
	for (i = 1; i < l; ++i) {
        if (i <= r)
            z[i] = min(r - i, z[i - j]);
        else
            z[i] = 0;
		for ( ; st[z[i]] == st[i + z[i]]; ++z[i]);
		if (i + z[i] > r) {
			r = i + z[i];
			j = i;
		}
	}
}
// Lyndon decomposition and minimal cyclical shift search in O(n)
string min_cyclic_shift(string s) {
	s += s;
	int n = (int) s.length();
	int i = 0, ans = 0;
	while (i < n / 2) {
		ans = i;
		int j = i + 1, k = i;
		while (j < n && s[k] <= s[j]) {
			if (s[k] < s[j]) k = i;
			else ++k;
			++j;
		}
		while (i <= k)  i += j - k;
	}
	return s.substr(ans, n / 2);
}
// Palindromes O(N): d1[i] - count of odd palindromes with center i, 
int l = 0, r = -1; //d2[i] - even palindromes with "right" center i
for (int i = 0; i < n; ++i) {
	int k = (i > r ? 1 : min(d1[l + r - i], r - i));
	while (i + k < n && i - k >= 0 && s[i + k] == s[i - k]) ++k;
	d1[i] = k--;
	if (i + k > r) {
		l = i - k;
		r = i + k;
	}
}
l = 0, r = -1;
for (int i = 0; i < n; ++i) {
	int k = (i > r ? 0 : min(d2[l + r - i + 1], r - i + 1)) + 1;
	while (i + k - 1 < n && i - k >= 0 && s[i + k - 1] == s[i - k]) ++k;
	d2[i] = --k;
	if (i + k - 1 > r) {
		l = i - k;
		r = i + k - 1;
	}
}