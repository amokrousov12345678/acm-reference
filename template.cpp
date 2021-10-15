#include <bits/stdc++.h>
//set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -fsanitize=address -fsanitize=undefined")
#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")//Comment optimisations for interactive problems (use endl)
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")
using namespace std;
#ifdef ILIKEGENTOO
void E(){}template<class A,class...B>void E(A $,B..._){cerr<<' '<<$;E(_...);}
# define E($...) E(#$,'=',$,'\n')
#else
# define E($...)
#endif
#define Sz(x) (int((x).size()))
#define All(x) begin(x),end(x)
using ll = long long;
const int infI = 0x3f3f3f3f;
const int infLL = 0x3f3f3f3f3f3f3f3f;

#include <ext/pb_ds/assoc_container.hpp> // Main file
#include <ext/pb_ds/tree_policy.hpp> // Contains tree_order_statistics_node_update
using namespace __gnu_pbds;
typedef tree<int, __gnu_pbds::null_type, less<int>, __gnu_pbds::rb_tree_tag, 
		__gnu_pbds::tree_order_statistics_node_update>
        ordered_set; //find_by_order and order_of_key
uset.max_load_factor(hashTableLoadRate); //0.25 for faster shesh tables

std::random_device rd; std::mt19937 rng{rd()}; std::mt19937 randMT(rng);//uniform int [a;b]
int blessRng(int a, int b) { std::uniform_int_distribution<int> rang(a, b); return rang(randMT);}
int main() {
	cin.tie(0); cout.tie(0);
    ios_base::sync_with_stdio(false);//AHTUNG: esli zabudesh, TL ottarabanit
	cout << setprecision(15) << fixed;
	//((float)(clock() - t0)) / CLOCKS_PER_SEC
    return 0;
}
//For tight ML: c++ io, vectors have very low impact (about 100kb). But compiler version/bitness is important factor
//DFS may cause ML (big recursion)
//euler formula: n - m + f = 2: m <= 3*n-6 - for planarity; n-verts, m - edges, f - faces (with outer)
//simpson integrate: (f(x0)+4f(x1)+2f(x2)+4f(x3)+..+4f(x_(2n-1))+f(x_(2n)))*h/3 points [x0..x_(2n)] evenly
//f(n) = sum(d|n) g(d) <=> g(n) = sum(d}n)(mu(d)*f(n/d)) (mobius inversion). n = sum(d|n)phi(d)
