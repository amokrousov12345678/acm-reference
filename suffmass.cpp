// ************** Suffix array construction (checked) *********************
char st[maxl]; //input str ('\0') terminated and with $ at end
int ind[maxl]; //res (sufs numbers in lex order)
int nind[maxl], cnt[maxl], cc, log_2, l;
int val[20][maxl]; //val[k][n] - id of equivalence class of [pos; pos+2^k) among strings of such len

void build_mass() {
	int i, j;
	l = strlen(st);
	for(log_2 = 0; (1 << log_2) < l; ++log_2) ;
	memset(cnt, 0, sizeof(cnt));
	memset(val, 0, sizeof(val));
	for (i = 0; i < l; ++i) cnt[st[i] + 1]++;
	for (i = 0; i < 128; ++i) cnt[i + 1] += cnt[i];
	for (i = 0; i < l; ++i) ind[cnt[st[i]]++] = i;
	for (i = 1; i < l; ++i) val[0][ind[i]] = val[0][ind[i-1]] + (st[ind[i]] != st[ind[i - 1]]);
	for (j = 0; j < log_2; ++j) { //use val[j % 2], val[(j+1)%2] if need to save memory
		cnt[0] = 0; cnt[cc = 1] = 1;
		for (i = 1; i < l; ++i)
			if (val[j][ind[i]] == val[j][ind[i - 1]]) ++cnt[cc];
			else cnt[++cc] = 1;
		for (i = 0; i < cc; ++i) cnt[i + 1] += cnt[i];
		int sz = (1 << j);
		for (i = 0; i < l; ++i) nind[ cnt[val[j][(ind[i] - sz + l) % l]]++ ] = (ind[i] - sz + l) % l;
		memcpy(ind, nind, sizeof(int) * l);
		for(i = 1; i < l; ++i) 
			val[j + 1][ind[i]] = val[j + 1][ind[i - 1]] + 
			(val[j][ind[i]] != val[j][ind[i - 1]] || val[j][(ind[i] + sz) % l] != val[j][(ind[i - 1] + sz) % l]);
	}
}

int lcp(int i, int j) { //suf start at i and at j
	int ans = 0, k;
	for(k = log_2 - 1; k >= 0; --k)
		if(val[k][i] == val[k][j]) {
			i = (i + (1 << k)) % l;
			j = (j + (1 << k)) % l;
			ans += (1 << k);
		}
	return ans;
}

int lcps[maxl], lcpLifts[20][maxl], pos[maxl]; //pos - inv of ind[], lcps - of adj in lex order
void build_lcp() {
	for (int i = 0; i < l; ++i) pos[ind[i]] = i;
	int k = 0;
	for (int i = 0; i < l; ++i) {
		if (k > 0) k--;
		if (pos[i] == l - 1) {
			lcps[l - 1] = -1; k = 0; continue;
		} else {
			int j = ind[pos[i] + 1];
			while (max(i + k,j + k) < l && st[(i + k) % l] == st[(j + k) % l]) ++k;
			lcps[pos[i]] = k;
		}
	}
    //just sparse table, which can be used for any other purposes (lcps[i] - ref array)
	for (int i=0;i<l;i++) lcpLifts[0][i] = lcps[i];
	int sz = 2;
	for (int i=1;i<=log_2;i++) {
		for (int j = 0; j + sz <= l; j++)
			lcpLifts[i][j] = min(lcpLifts[i - 1][j], lcpLifts[i - 1][j + sz / 2]);
		sz <<= 1;
	}
}

int lcpKasai(int i, int j) { //lcp.rmq[pos(i);pos(j)) O(1) sparse table
	i = pos[i]; j = pos[j]; if (i>j) swap(i,j);//don't forget SWAP
	if (i==j) return l;
	int t = __lg(j-i);
	return min(lcpLifts[t][i], lcpLifts[t][j - (1 << t)]);
}