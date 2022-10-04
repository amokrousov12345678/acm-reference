// ********************* Match in bipartite graph **********************
vector<int> g[maxn];//MUST BE undirected if edge (a,b) then should be (b,a)
int used[maxn];
int match[maxn];//match vert id (correct for verts from both parts)
int cc; //num of found edges
//|MinVertexCover| = N - |MaxIndSet| (complementary to any vert cover is indepedent set)
//|MinVertexCover| = |Matching| (for bipartite).
//How to build: matching edges as R->L, other L->R, dfs from all FREE vertices in L
//+,- about visited. Then L- union R+ is ans (easy to prove from L+R- is ind.set)
//Weighted MinVC - min cut(=maxflow) with cap to verts=their cost and INF cap for edges
//If weight negative-add weight to ans, but consider cap as 0
bool dfs(int v) {
    if (used[v] == cc) return false;
    used[v] = cc;
    for (auto& it: g[v]) {
        if (match[it]==-1 || dfs(match[it])) {
            match[v] = it; match[it] = v; return true;
        }
    }
    return false;
}
int kuhn(int n) {//try to match [1..n] (enough only one part, but may try ALL)
    memset(used, -1, n*sizeof(int));//MUST clear for full size of G, not of part
    memset(match, -1, n*sizeof(int));//SAME
    cc = 0;
    for (int i=0;i<n;i++) {//may choose in other order (for example based on weight)
        if (match[i]!=-1) continue;
        if (dfs(i)) cc++;
    }
    return cc;
}