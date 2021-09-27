// **************************** Range trees ****************************
const int inf = 0x3f3f3f3f;
typedef int val_t;
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
//if query op is addition, we need store len[v] to scale impact of stored addition
struct rmq_rm_t { // ***** RMQ with interval modification support *****
	val_t ans[size << 1], add[size << 1];
	rmq_rm_t() { }
	rmq_rm_t(val_t *a) {
		memset(add, 0, sizeof(add)); // ! you may need something else here
		memcpy(ans + size, a, sizeof(val_t) * size);
		for (int i = size - 1; i > 0; --i) ans[i] = min(ans[i << 1], ans[(i << 1) + 1]);
	}
	void modify(int l, int r, val_t val) {
		if (l > r) return;
		l += size; r += size;
		add[l] += val; if (l < r) add[r] += val;
		while (l > 0) {
			if ((l & 1) == 0 && (l + 1) < r) add[l + 1] += val;
			if ((r & 1) == 1 && (r - 1) > l) add[r - 1] += val;
			l >>= 1; r >>= 1;
			ans[l] = min(ans[l << 1] + add[l << 1], ans[(l << 1) + 1] + add[(l << 1) + 1]);
			ans[r] = min(ans[r << 1] + add[r << 1], ans[(r << 1) + 1] + add[(r << 1) + 1]);
		}
	}
	val_t query(int l, int r) {
		if (l > r) return inf;
		l += size; r += size;
		val_t lans = ans[l] + add[l], rans = ans[r] + add[r];
		//ll left = len[l], right = l<r ? len[r] : 0;
		while (l > 0) {
			if ((l & 1) == 0 && (l + 1) < r) lans = min(lans, ans[l + 1] + add[l + 1]);
			//if ((l & 1) == 0 && (l + 1) < r) left += len[l+1];
			if ((r & 1) == 1 && (r - 1) > l) rans = min(rans, ans[r - 1] + add[r - 1]);
			//if ((r & 1) == 1 && (r - 1) > l) right += len[r-1];
			l >>= 1; r >>= 1;
			lans += add[l]; rans += add[r]; //lans += add[l] * left; rans += add[r] * right;
		}
		return min(lans, rans);
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
  
    val_t applyUpdate(val_t val, val_t push) {};
    val_t rqOp(val_t lhs, val_t rhs) {}
    val_t combineUpdates(val_t cur, val_t fromUp) {return cur+fromUp;};

    void build(int* vals, int n, int v = 1, int l = 0, int r = size-1) {
        if (l==r) {
            mas[v] = l<n ? vals[l] : neutral;
            def[v] = neutralUpdate;
            return;
        }
        int m = (l+r)/2;
        build(vals, n, 2*v, l, m);
        build(vals, n, 2*v+1, m+1, r);
        mas[v] = rqOp(applyUpdate(mas[2*v], def[2*v]), applyUpdate(mas[2*v+1], def[2*v+1]));
        def[v] = neutralUpdate;
    }
    void push(int v) { /*push(Node*&)*/
        def[2*v] = combineUpdates(def[2*v], def[v]);
        def[2*v+1] = combineUpdates(def[2*v+1], def[v]);
        def[v] = neutralUpdate;
    }

    void update(int rq_l, int rq_r, val_t val, int v = 1, int l = 0, int r = size-1) {
        if (l > rq_r || r < rq_l) return;
        if (rq_l <= l && r <= rq_r) {
            def[v] = combineUpdates(def[v], val);
            return;
        }
        push(v);
        int m = (l+r)/2;
        update(rq_l, rq_r, val, 2*v, l, m);
        update(rq_l, rq_r, val, 2*v+1, m+1, r);
        mas[v] = rqOp(applyUpdate(mas[2*v], def[2*v]), applyUpdate(mas[2*v+1], def[2*v+1]));
    }

    /*pair<val_t, Node*>*/val_t query(int rq_l, int rq_r, int v = 1, int l = 0, int r = size-1) {
        if (l > rq_r || r < rq_l) return neutral;
        if (rq_l <= l && r <= rq_r) {
            return applyUpdate(mas[v], def[v]);
        }
        push(v);
        int m = (l+r)/2;
        auto ans = rqOp(query(rq_l, rq_r, 2*v, l, m), query(rq_l, rq_r, 2*v+1, m+1, r));
        mas[v] = rqOp(applyUpdate(mas[2*v], def[2*v]), applyUpdate(mas[2*v+1], def[2*v+1]));
        return ans;
    }
};