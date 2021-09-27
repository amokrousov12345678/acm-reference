// ************** Construction of suffix automata in O(n) **************
const int maxn = 2e5; //states count
const int maxa = 53; //alphabet size
char s[maxl]; //input string ('\0' terminated)

//TODO: add firstPos
struct node { //vertex correspond to class of string with equal endpos sets
	int next[maxa], suff, len, end;//end - if terminal, len - longest string in endpos equality class 
	//suff - link to vertex of class of longest suffix which still isn't in same class as string
	//to find all occurences of string, do dfs on reverse suf links tree (+don't print clones, but traverse them)
	int firstPos;//first endpos of vertex
} nodes[maxn];
int root, cnt;
inline int new_node(int len) {
	memset(nodes + cnt, -1, sizeof(node));
	nodes[cnt].len = len;
	nodes[cnt].end = 0;
	return cnt++;
}
void extend(char c, int &last) {
	int nlast = new_node(nodes[last].len + 1), p;
	nodes[nlast].firstPos = nodes[nlast].len - 1;
	for (p = last; p >= 0 && nodes[p].next[numc(c)] == -1; p = nodes[p].suff)
		nodes[p].next[numc(c)] = nlast;
	last = nlast;
	if (p < 0) {
		nodes[nlast].suff = root;
		return;
	}
	int q = nodes[p].next[numc(c)];
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
	int l = strlen(s), i, last;
	last = root = new_node(0);
	for (i = 0; i < l; ++i) extend(s[i], last);
	for (i = last; i >= 0; i = nodes[i].suff) nodes[i].end = 1;
}