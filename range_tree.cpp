// **************************** Range trees ****************************
const int size = 1 << 17;//MUST BE power of 2
struct rmq_t { // ********* RMQ without interval modification *********
	val_t mas[size << 1];
	rmq_t() {}
	rmq_t(val_t *a) {
		memcpy(mas + size, a, sizeof(val_t) * size);
		for (int i = size - 1; i > 0; --i) mas[i] = min(mas[i << 1], mas[(i << 1) + 1]);
	}
	void modify(int ind, val_t val) {
		ind += size;
		mas[ind] += val;
		for (ind >>= 1; ind > 0; ind >>= 1) mas[ind] = min(mas[ind << 1], mas[(ind << 1) + 1]);
	}
	val_t query(int l, int r) {
		if (l > r) return inf;
		l += size; r += size;
		val_t ans = mas[l]; if (l<r) ans = min(ans, mas[r]);
		for ( ; l < r; l >>= 1, r >>= 1) {
			if ((l & 1) == 0 && (l + 1) < r) ans = min(ans, mas[l + 1]);
			if ((r & 1) == 1 && (r - 1) > l) ans = min(ans, mas[r - 1]);
		}
		return ans;
	}
};
// **************************** Fenvik tree ****************************
struct fenvic_t { //for multiple dimensions loops by each index
    val_t mas[size];
    fenvic_t() {};
    void modify(int ind_, val_t val) {
        for (int ind = ind_; ind < size; ind = (ind | (ind+1))) mas[ind] += val;
    }
    val_t query(int ind_) {
        if (ind_<0) return val_t();
        val_t ans = val_t();
        for (int ind = ind_; ind >= 0; ind = (ind & (ind+1)) - 1) ans += mas[ind];
        return ans;
    }
    val_t query(int l, int r) {
        return query(r) - query(l-1);
    }
};
//for recursive, we need recreate vertices on ANY change, even when push (otherwise children are desync)
struct segtree {
    val_t mas[4*size]; const static val_t neutral = 0;
    pval_t def[4*size]; const static pval_t neutralUpdate = 0;
    val_t applyUpdate(val_t val, pval_t push) {};//apply update to  vertex
    val_t rqOp(val_t lhs, val_t rhs) {}//op for range query
    pval_t combineUpdates(pval_t cur, pval_t fromUp) {return cur+fromUp;};//apply def[par] into def[v]

    void build(int* vals, int n, int v = 1, int l = 0, int r = size-1) {
        if (l==r) {mas[v] = l<n ? vals[l] : neutral; def[v] = neutralUpdate; return;}
        int m = (l+r)/2;
        build(vals, n, 2*v, l, m); build(vals, n, 2*v+1, m+1, r);
        mas[v] = rqOp(applyUpdate(mas[2*v], def[2*v]), applyUpdate(mas[2*v+1], def[2*v+1]));
        def[v] = neutralUpdate;
    }
    void push(int v) { /*push(Node*&)*/
        def[2*v] = combineUpdates(def[2*v], def[v]);
        def[2*v+1] = combineUpdates(def[2*v+1], def[v]); def[v] = neutralUpdate;
    }
    void update(int rq_l, int rq_r, pval_t val, int v = 1, int l = 0, int r = size-1) {
        if (l > rq_r || r < rq_l) return;
        if (rq_l <= l && r <= rq_r) {def[v] = combineUpdates(def[v], val); return;}
        push(v); int m = (l+r)/2;
        update(rq_l, rq_r, val, 2*v, l, m); update(rq_l, rq_r, val, 2*v+1, m+1, r);
        mas[v] = rqOp(applyUpdate(mas[2*v], def[2*v]), applyUpdate(mas[2*v+1], def[2*v+1]));
    }
    /*pair<val_t, Node*>*/val_t query(int rq_l, int rq_r, int v = 1, int l = 0, int r = size-1) {
        if (l > rq_r || r < rq_l) return neutral;
        if (rq_l <= l && r <= rq_r) return applyUpdate(mas[v], def[v]);
        push(v); int m = (l+r)/2;
        auto ans = rqOp(query(rq_l, rq_r, 2*v, l, m), query(rq_l, rq_r, 2*v+1, m+1, r));
        mas[v] = rqOp(applyUpdate(mas[2*v], def[2*v]), applyUpdate(mas[2*v+1], def[2*v+1]));
        return ans;
    }
    //TODO: add bin search (from left)
};
//LiChao tree: store line in vert. On process vert: choose best (min/max) among newF and mas[v] in pnt m
//Based on comparison in pnt l, in one half newF never could be optimal, descend to other
//Leafs trivially updated. Answer to query: best in given point on path to leaf
struct segtree {//persistent (WIHTOUT group ops, for them adopt upper version (THINK 150 TIMES!!!))
    static val_t rqOp(val_t lhs, val_t rhs) {return lhs+rhs;}
    struct Node {
        Node *l, *r; val_t val; Node(): l(nullptr), r(nullptr), val(neutral) {};
        Node(Node *l, Node *r): l(l), r(r) { val = rqOp((l ? l->val : neutral), (r? r->val : neutral));}
    };
	
    Node* build(int l = 0, int r = size-1) {
        if (l==r) return new Node();
        int m = (l+r)/2; return new Node(build(l, m), build(m+1, r));
    }
    Node* modify(int pos, val_t newVal, Node* v, int l  = 0, int r = size-1) {
        if (l==r) {auto tmp = new Node(*v); tmp->val = newVal; return tmp;}
        int m = (l+r)/2;
        auto lCh = pos<=m ? modify(pos, newVal, v->l, l, m) : v->l;
        auto rCh = pos>m ? modify(pos, newVal, v->r, m+1, r) : v->r;
        return new Node(lCh, rCh);
    }
    val_t query(int rq_l, int rq_r, Node* v, int l = 0, int r = size-1) {
        if (l > rq_r || r < rq_l) return neutral;
        if (rq_l <=l && r<=rq_r) return v->val;
        int m = (l+r)/2;
        val_t res = rqOp(query(rq_l, rq_r, v->l, l, m), query(rq_l, rq_r, v->r, m+1, r));
        return res;
    }
};//for lazy on big range create children for vertex on demand