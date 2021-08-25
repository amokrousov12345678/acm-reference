// **************************** Range trees ****************************
//GRABLES:
//if op isn't idempotent, must check l!=r before adding mas[r]
//if query op is addition, we need store len[v] to scale impact of stored addition

const int inf = 0x3f3f3f3f;
typedef int val_t;
const int size = 1 << 17;

//for recursion array mas[size<<2] and ALL methods should be recursive (modify, etc, etc)
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
		val_t ans = min(mas[l], mas[r]);
		for ( ; l < r; l >>= 1, r >>= 1) {
			if ((l & 1) == 0 && (l + 1) < r) ans = min(ans, mas[l + 1]);
			if ((r & 1) == 1 && (r - 1) > l) ans = min(ans, mas[r - 1]);
		}
		return ans;
	}
};

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
		add[l] += val;
		if (l < r) add[r] += val;
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
