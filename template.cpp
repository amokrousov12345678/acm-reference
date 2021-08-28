#include <bits/stdc++.h>

#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")//Comment optimisations for interactive problems (use endl)
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

#include <ext/pb_ds/assoc_container.hpp> // Общий файл.
#include <ext/pb_ds/tree_policy.hpp> // Содержит класс tree_order_statistics_node_update
using namespace __gnu_pbds;
typedef tree<
        int,
        __gnu_pbds::null_type,
        less<int>,
        __gnu_pbds::rb_tree_tag,
        __gnu_pbds::tree_order_statistics_node_update>
        ordered_set; //find_by_order and order_of_key

using namespace std;
#ifdef ILIKEGENTOO
void E(){}template<class A,class...B>void E(A $,B..._){cerr<<' '<<$;E(_...);}
# define E($...) E(#$,'=',$,'\n')
#else
# define E($...)
#endif
#define Sz(x) (int((x).size()))
#define All(x) begin(x),end(x)

typedef double flt;
typedef long long int64;

const int infI = 0x3f3f3f3f;
const int infLL = 0x3f3f3f3f3f3f3f3f;

int main() {
	cin.tie(0);
    cout.tie(0);
    ios_base::sync_with_stdio(false);
	cout << setprecision(15) << fixed;
	//((float)(clock() - t0)) / CLOCKS_PER_SEC
    return 0;
}
