#include <bits/stdc++.h>
//set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -fsanitize=address -fsanitize=undefined")
//#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math")
//#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,fma")//drop avx and fma if RE#1
//(on Yandex, NSUts remove also avx)
////#pragma GCC optimize("O3,no-stack-protector,unroll-loops") FOR GEOM INSTEAD OF 1st
////#pragma comment(linker, "/STACK:36777216") (for Studio if get ML on GCC)
using namespace std;
#ifdef ILIKEGENTOO
void E(){}template<class A,class...B>void E(A $,B..._){cerr<<' '<<$;E(_...);}
# define E($...) E(#$,'=',$,'\n')
#else
# define E($...)
#endif
#define Sz(x) (int((x).size()))
#define assertTL(x) {if (!(x)) while(1);};
using ll = long long;
const int infI = 0x3f3f3f3f;
const ll infLL = 0x3f3f3f3f3f3f3f3f;

struct hash_pair { //to allow unordered_map<pair<int, int>, bool, hash_pair>
	template <class T1, class T2> size_t operator()(const pair<T1, T2>& p) const {
		auto hash1 = hash<T1>{}(p.first); auto hash2 = hash<T2>{}(p.second); return hash1 ^ hash2;}};
// __builtin_popcount

std::random_device rd; std::mt19937 rng{rd()}; std::mt19937 randMT(rng);//uniform int [a;b]
int blessRng(int a, int b) { std::uniform_int_distribution<int> rang(a, b); return rang(randMT);}
int main() {    
	cin.tie(0); cout.tie(0);
    ios_base::sync_with_stdio(false);//WARNING: esli zabudesh, TL (krodet'sya)
	cout << setprecision(15) << fixed;
	//((float)(clock() - t0)) / CLOCKS_PER_SEC
    return 0;
}
const int MAX_MEM = 1e8; int mpos = 0; //Fast allocator (danger, use ONLY as last chance!)
char mem[MAX_MEM]; inline void * operator new ( size_t n ) {char *res = mem + mpos; mpos += n;
assert(mpos <= MAX_MEM); return (void *)res; }; inline void operator delete ( void * ) { }
//For tight ML: c++ io, vectors have very low impact (about 100kb). But COMPILER VERSION/BITNESS
//TL: vector<vector<>> makes A LOT heap allocations, so may be slower 10 times than just vector<int>
//Or use vector<int> arr[maxn] and cleanup. Try to compress coords to avoid use of map
//is important factor. DFS may cause ML (big recursion). On rect grid use BFS
//Euler formula: n - m + f = 1+cc: m <= 3*n-6 - for planarity; n-verts, m - edges, f - faces (with outer)
//simpson integrate: (f(x0)+4f(x1)+2f(x2)+4f(x3)+..+4f(x_(2n-1))+f(x_(2n)))*h/3 points [x0..x_(2n)] evenly
//f(n) = sum(d|n) g(d) <=> g(n) = sum(d}n)(mu(d)*f(n/d)) (mobius inversion). n = sum(d|n)phi(d)
//Kotlovan numbers: C(n) = sum(i=0;i<=n-1)(C(i)*C(n-1-i)) = 1/(n+1)*Binom(2n,n). 
//Number of correct braces sequence len 2n
//Newton method: x0, x_(i+1) = x_i - f(x_i)/f'(x_i). Need good start approx. Can used to calc sqrt
//Stable matching: each vert has prefer list order, match to avoid conflict:
//(A B) and (a b) matched, but A prefers b more than a, b prefers A more than B
//Solution: free vert of left part tries best in his opinion not checked option. 
//Second part accept rebind if new option better. Continue until all matched or pref lists exhausted.
//Pik formula: S = I+B/2 -1. S-area, I-count int points inside, B-count int points of bound. S must be >0
//std::string_stream is cool thing to tokenize string)
//L1 (|x|+|y|) metric with replace x'=x+y, y'=x-y transformed to Linf (max(|x|,|y|))
System.setIn(new FileInputStream("input.txt"));//LONG ARITHM, BigInteger immutable)
Scanner scanner = new Scanner(new BufferedReader(new InputStreamReader(System.in)));

#include <ext/rope> //alternative for treap from STL (works on GCC)
using namespace __gnu_cxx; rope<int> v(n, 0); //rope<int> v(n) DOESN'T WORK AS EXPECTED
//.substr(l, len) .insert(l, seg) .erase(l, len) .mutable_reference_at(pos) .mutable_begin .mutable_end

#include <ext/pb_ds/assoc_container.hpp> // Main file for ordered set
#include <ext/pb_ds/tree_policy.hpp> // Contains tree_order_statistics_node_update
using namespace __gnu_pbds;
typedef tree<int, __gnu_pbds::null_type, less<int>, __gnu_pbds::rb_tree_tag, 
		__gnu_pbds::tree_order_statistics_node_update>
        ordered_set; //find_by_order and order_of_key (0-based, don't work as multiset)
uset.max_load_factor(hashTableLoadRate); //0.25 for faster shesh tables