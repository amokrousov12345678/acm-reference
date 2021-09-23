#define numc(x) ((int)((unsigned char)x)) //if alphabet contains negative codes
const int maxa = 53; //alphabet capacity
const int maxn = 2e5; //verts count
const int maxw = 2e5; //word count
struct node {
	int end, suff, cSuff, next[maxa]; //end - id of one of word in vertex
	//cSuff - compressed suf links (if we need poses of all occurences)
} nodes[maxn];
int cnt = 0, root = 0;
vector<int> same_words[maxw]; //indexes of words equal to given (written as node[x].end)
int word_lens[maxw];

int new_node() {
	memset(nodes + cnt, -1, sizeof(node));
	return cnt++;
}
void init() { new_node(); } //MUST be called before use of trie
void add_word(const char *st, int num) {
	assert(cnt>0); word_lens[num] = strlen(st); int v; 
	for (v = root; *st; ++st) {
		if (nodes[v].next[numc(*st)] == -1) nodes[v].next[numc(*st)] = new_node();
		v = nodes[v].next[numc(*st)];
	}
	if (nodes[v].end == -1) nodes[v].end = num;
	same_words[nodes[v].end].push_back(num);
}
int q[maxn]; //queue for bfs through vertices
void aho_corasik() { //(only suffix links) automaton construction
	assert(cnt>0); 
	nodes[root].suff = root;
	int l, r = 0, i, u, nx;
	for (i = 0; i < maxa; ++i)
		if ((nx = nodes[root].next[i]) != -1) {
			nodes[nx].suff = root;
			q[r++] = nx;
		} else nodes[root].next[i] = root;
	for (l = 0; l < r; ++l) {
		for (i = 0; i < maxa; ++i) if ((nx = nodes[q[l]].next[i]) != -1) {
			for (u = nodes[q[l]].suff; nodes[u].next[i] == -1; u = nodes[u].suff);
			nodes[nx].suff = nodes[u].next[i];
			q[r++] = nx;
		}
	}	
}
int was[maxn], pres[maxw];//pres - if given word found
void search(const char *st) {
	memset(was, 0, sizeof(int) * cnt);
	for (int v = root; *st; ++st) {
		for ( ; nodes[v].next[numc(*st)] == -1; v = nodes[v].suff) ;
		v = nodes[v].next[numc(*st)];
		for (int u = v; u != root; u = nodes[u].suff)
			if (was[u]++) break;
			else if (nodes[v].end != -1)
				for (int i = 0; i < same_words[nodes[v].end].size(); ++i)
					pres[same_words[nodes[v].end][i]] = 1;
	}
}
int getShortSuff(int v) {//lazy eval of cSuff (for traverse occurences ends in each pos)
	if (nodes[v].cSuff == -1) {
		if (nodes[nodes[v].suff].end != -1) nodes[v].cSuff = nodes[v].suff;
		else if (nodes[v].suff == root) nodes[v].cSuff = root; 
		else nodes[v].cSuff = getShortSuff(nodes[v].suff);
	}
	return nodes[v].cSuff;
}
int v = root;
for (int j=0;j<Sz(s);j++) {
	for ( ; nodes[v].next[numc(s[j])] == -1; v = nodes[v].suff);
	v = nodes[v].next[numc(s[j])];
	int u = nodes[v].end != -1 ? v : getShortSuff(v);
	for (; u != root; u = getShortSuff(u)) {
		OCCURENCE(j+1-word_lens[nodes[u].end]);
	}
}