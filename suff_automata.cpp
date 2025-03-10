// ************** Construction of suffix automata in O(n) **************
const int maxn = 2e5; //max states count (at least 2*maxl)
const int maxa = 53; //alphabet size
char s[maxl]; //input string ('\0' terminated)

struct node { //vertex correspond to class of string with equal endpos sets
	int next[maxa], suff, len, end;//end - if terminal, len - longest string in endpos equality class 
	//suff - link to vertex of class of longest suffix which not in same class as string
    //strings is class is suffixes of longest string with len in range (len(suf(v));len(v)]
	//to find all occurences of string, do dfs on reverse suf links tree (+DON'T PRINT clones, but TRAVERSE them)
	int firstPos;//first endpos of vertex
} nodes[maxn];
int root, cnt;
inline int new_node(int len) {
	memset(nodes + cnt, -1, sizeof(node));
	nodes[cnt].len = len; nodes[cnt].end = 0; nodes[cnt].firstPos = 0;
	return cnt++;
}
void extend(char c, int &last) {
	int nlast = new_node(nodes[last].len + 1), p;
	nodes[nlast].firstPos = nodes[nlast].len - 1;
	for (p = last; p >= 0 && nodes[p].next[numc(c)] == -1; p = nodes[p].suff)
		nodes[p].next[numc(c)] = nlast;
	last = nlast;
	if (p < 0) {//len(max existing suf of s+c) = 0
		nodes[nlast].suff = root;
		return;
	}
	int q = nodes[p].next[numc(c)]; //len(max existing suf of s+c)=nodes[p].len+1
	if (nodes[q].len == nodes[p].len + 1) {
		nodes[nlast].suff = q;
		return;
	}
	int nq = new_node(nodes[p].len + 1); //part which creates clone
	memcpy(nodes[nq].next, nodes[q].next, sizeof(nodes[q].next));
	nodes[nlast].suff = nq;
	nodes[nq].suff = nodes[q].suff;
	nodes[nq].firstPos = nodes[q].firstPos;
	nodes[q].suff = nq;
	for ( ; p >= 0 && nodes[p].next[numc(c)] == q; p = nodes[p].suff)
		nodes[p].next[numc(c)] = nq;
}
void create_automata() {
	int l = strlen(s), i, last; cnt = 0;
	last = root = new_node(0);
	for (i = 0; i < l; ++i) extend(s[i], last);
	for (i = last; i >= 0; i = nodes[i].suff) nodes[i].end = 1;
}//to build SA on several strings, connect them through "$" and don't go by "$"
string lcs (string s, string t) {//works only with root=0
	int v = 0,  l = 0, best = 0, bestpos = 0; //you MUST have built automaton on s
	for (int i=0; i<(int)t.length(); ++i) {
        while (v && (nodes[v].next[numc(t[i])]==-1)) {v = nodes[v].suff; l = nodes[v].len;}
        if (nodes[v].next[numc(t[i])]!=-1) {v = nodes[v].next[numc(t[i])]; ++l;}
        if (l > best) {best = l;  bestpos = i;}//now "l" - max len of lcs ended in t[i]
    } return t.substr (bestpos-best+1, best);
}