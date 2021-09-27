#include <bits/stdc++.h>

//set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -fsanitize=address -fsanitize=undefined")

#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")//Comment optimisations for interactive problems (use endl)
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

#include <ext/pb_ds/assoc_container.hpp> // Main file
#include <ext/pb_ds/tree_policy.hpp> // Contains tree_order_statistics_node_update
using namespace __gnu_pbds;
typedef tree<int, __gnu_pbds::null_type, less<int>, __gnu_pbds::rb_tree_tag, 
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

uset.max_load_factor(hashTableLoadRate); //0.25 for faster shesh tables

int main() {
	cin.tie(0);
    cout.tie(0);
    ios_base::sync_with_stdio(false);
	cout << setprecision(15) << fixed;
	//((float)(clock() - t0)) / CLOCKS_PER_SEC
    return 0;
}

//euler formula: n - m + f = 2, m <= 3*n-6
//simpson integrate: (f(x0)+4f(x1)+2f(x2)+4f(x3)+..+4f(x(2n-1))+f(x(2n)))*h/3
//For tight ML: c++ io, vectors have very low impact (about 100kb). But compiler version/bitness is important factor