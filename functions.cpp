// #include <bits/stdc++.h>
#include <iostream>
#include <map>
#include <cmath>
#include <ios>
#include <iomanip>
#include <algorithm>
#include <queue>
#include <stack>
#include <numeric>
// #include <windows.h>
using namespace std;
template<class T> using V = vector<T>;
template<class T> using VV = V<V<T>>;
template<class T> using VVV = V<VV<T>>;
template<class T1, class T2> using P = pair<T1, T2>;
using I = int;
using D = double;
using B = bool;
using C = char;
using S = string;
using LL = long long;
using LD = long double;
using ULL = unsigned long long;
using PII = P<I, I>;
using VPII = V<PII>;
using PLL = P<LL, LL>;
using VPLL = V<PLL>;
using VI = V<I>;
using VVI = VV<I>;
using VLL = V<LL>;
using VVLL = VV<LL>;
using VC = V<C>;
using VVC = VV<C>;
using VS = V<S>;
using VVS = VV<S>;
using VB = V<B>;
using VVB = VV<B>;
#define REP(type, i, n) for (type i = 0; i < (type)(n); ++i)
#define REP2(type, i, m, n) for (type i = (m); i < (type)(n); ++i)
#define REPR(type, i, n) for (type i = (n)-1; i >= 0; --i)
#define REPR2(type, i, m, n) for (type i = (n)-1; i >= (m); --i)
#define REPx(x, a) for(auto x : a)
#define ALL(a) a.begin(), a.end()
#define SORT(a) sort(ALL(a))
#define SORTR(a, type) sort(ALL(a), greater<type>())
#define SORTF(a, comp) sort(ALL(a), comp)
#define REVERSE(a) reverse(ALL(a))
#define SIZE(a, type) ((type)(a).size())
#define bit_search(bit, n) REP(LL, bit, 1LL<<(n))
#define bit_check(bit, i) ((bit>>(i)) & 1)
#define setpre(n) fixed << setprecision((n))
#define UNIQUE(a) do {SORT(a); (a).erase(unique(ALL(a)), (a).end());} while(0)
#define MAX(a) *max_element(ALL(a))
#define MIN(a) *min_element(ALL(a))
#define bisect_left(a, x) (lower_bound(ALL(a), (x)) - a.begin())
#define bisect_right(a, x) (upper_bound(ALL(a), (x)) - a.begin())
#define INPUT(a) REPx(&x, a) cin >> x;
#define INPUT2(a) REPx(&x, a) INPUT(x);
#define INPUTP(a) REPx(&x, a) cin >> x.first >> x.second;
#define ENDL cout << endl;

const int INF = 2e9;
// const LL INF = 9e18;
const LL MOD = 1e9+7;

template<class T> using PRIORITY_QUEUE = priority_queue< T, V<T>, greater<T> >;
template<class T> inline bool chmin(T &a, T b){if (a > b) {a = b; return true;} return false;}
template<class T> inline bool chmax(T &a, T b){if (a < b) {a = b; return true;} return false;}
template<class T> inline void debug1(V<T> A){REP(int, i, SIZE(A, int)){if (A[i] == INF) cout << "I ";else cout << A[i] << " ";}ENDL}
template<class T> inline void debug2(VV<T> A){REP(int, i, SIZE(A, int)){REP(int, j, SIZE(A[i], int)){if (A[i][j] == INF) cout << "I "; else cout << A[i][j] << " ";}ENDL}}
/*
ファイルから入力を受け付ける。(cygwin では ./aaa < test.text のようにコマンドを使えばいい)
  #include <fstream>
  ifstream in("test.txt");
  cin.rdbuf(in.rdbuf());
*/


template<class T>
T GCD(T a, T b)
{
  while(a % b)
  {
    T r = a % b;
    a = b;
    b = r;
  }
  return b;
}

template<class T>
T LCM(T a, T b)
{
  return (a * b) / GCD(a, b);
}

template<class T>
T gcd(T a, T b)
{
  // GCD 再帰関数バージョン
  if (a % b == 0) return b;
  return gcd(b, a % b);
}

template<class T>
T Pow(T a, T b)
{
  // a^b を計算する。
  T res = 1;
  while (b > 0)
  {
    if (b & 1) res = res * a;
    a = a * a;
    b >>= 1;
  }
  return res;
}

template<class T>
T comb(T n, T k)
{
  /*
  nCkを計算する。
            n!        n(n-1)...(n-k+1)     n     n-1     n-2           n-k+1
  nCk = ---------- = ------------------ = --- * ----- * ----- * ... * -------
         k!(n-k)!        k(k-1)...1        1      2       3              k
  と分解すると、右辺の第k項目までの積は整数になるので、逐次割り算計算できてオーバフローになりにくくなる。
  (∵連続するk整数の積はk!で割り切れる)
  */
  if (n < k) return 0;
  if (k > n/2) k = n - k;
  T res = 1;
  REP(T, i, k) res = res * (n-i) / (i+1);
  return res;
}

template<class T>
void make_nCk(VV<T>& indexes, T n, T k)
{
  /*
  nCkを列挙する。VVI indexes に a[i] = i(i=1~n) という数列の組合せを格納する。(n種類の数字からk個選ぶ組み合わせ)
  indexes のサイズをあらかじめ決定するために、nCkを計算する関数 int comb() を宣言しておき、
  VVI indexes(comb(n, k), VI(k));
  と宣言する。
  */
  T idx = 0, ct = 0;
  bit_search(bit, n)
  {
    ct = 0;
    REP(T, i, n) if (bit_check(bit, i)) ct++;
    if (ct != k) continue;
    ct = 0;
    REP(T, i, n) if (bit_check(bit, i)) indexes[idx][ct++] = i+1;
    if (++idx == SIZE(indexes, T)) return;
  }
}

template<class T>
void make_nHk(VV<T>& indexes, const T n, const T k, T s, T depth)
{
  /*
  nHkを列挙する。VVI indexes に a[i] = i(i=1~n) という数列の重複組合せを格納する。(n種類の数字から重複を許してk個選ぶ組み合わせ)
  indexes のサイズをあらかじめ決定するために、nCkを計算する関数 int comb() を宣言しておき、
  VVI indexes(comb(k+n-1, k), VI(k));
  と宣言する。(nHk = k+n-1Ck)
  引数の s, depth にははじめ、0を入れる。
  */
  static T idx = 0;
  static V<T> A(k);
  if (depth == k)
  {
    REP(T, i, k) indexes[idx][i] = A[i];
    idx++;
    return;
  }
  REP2(T, i, s, n)
  {
    A[depth] = i+1;
    make_nHk(indexes, n, k, i, depth+1);
  }
}

// nCkの列挙
template<class T>
struct Combination
{
  VV<T> indexes;

  Combination(const T& n, const T& k, const bool& is_sorted=false): indexes(comb(n, k), V<T>(k))
  {
    make_nCk(n, k);
    if (is_sorted) SORT(indexes); // is_sortedがtrueなら配列をソートする
  }

  T comb(T n, T k)
  {
    if (n < k) return 0;
    if (k > n/2) k = n - k;
    T res = 1;
    REP(T, i, k) res = res * (n-i) / (i+1);
    return res;
  }

  void make_nCk(const T& n, const T& k)
  {
    T idx = 0, ct = 0;
    bit_search(bit, n)
    {
      ct = 0;
      REP(T, i, n) if (bit_check(bit, i)) ct++;
      if (ct != k) continue;
      ct = 0;
      REP(T, i, n) if (bit_check(bit, i)) indexes[idx][ct++] = i+1;
      if (++idx == SIZE(indexes, T)) return;
    }
  }
};

// nHkの列挙
template<class T>
struct Homogeneous
{
  VV<T> indexes;

  Homogeneous(const T& n, const T& k): indexes(comb(k+n-1, k), V<T>(k))
  {
    make_nHk(n, k, 0, 0);
  }

  T comb(T n, T k)
  {
    if (n < k) return 0;
    if (k > n/2) k = n - k;
    T res = 1;
    REP(T, i, k) res = res * (n-i) / (i+1);
    return res;
  }

  void make_nHk(const T& n, const T& k, const T s, const T depth)
  {
    static T idx = 0;
    static V<T> A(k);
    if (depth == k)
    {
      REP(T, i, k) indexes[idx][i] = A[i];
      idx++;
      return;
    }
    REP2(T, i, s, n)
    {
      A[depth] = i+1;
      make_nHk(n, k, i, depth+1);
    }
  }
};
/*
int main(){
  int n, k;
  cin >> n >> k;
  Homogeneous<int> ho(n, k);
  REP(int, i, SIZE(ho.indexes, int)){
    REP(int, j, k) cout << ho.indexes[i][j] << " ";
    ENDL
  }
  return 0;
}
*/

// 順列の列挙(n=10が限界, 400~440 ms)
template<class T>
struct Permutation
{
  VV<T> indexes;

  Permutation(const T& n, const T& s=1): indexes(factorial(n), V<T>(n))
  {
    // s は最小の数字
    T idx = 0;
    V<T> v(n);
    iota(ALL(v), s);
    do
    {
      indexes[idx++] = v;
    } while (next_permutation(ALL(v)));
  }

  T factorial(const T& a)
  {
    T res = 1;
    REP2(T, i, 2, a+1) res *= i;
    return res;
  }
};

///////////////////////////////////////////////////////////////////////////////
// mod の計算用関数

template<class T>
T modpow(T a, T b, T mod)
{
  // a^b を mod で割った余りを計算する。
    T res = 1;
    while (b > 0)
    {
        if (b & 1) res = res * a % mod;
        a = a * a % mod;
        b >>= 1;
    }
    return res;
}

template<class T>
T modcomb(T n, T k, T mod)
{
  // nCkのmodを計算する。1<=n<=10^9、 1<=k<=10^7の場合に有効。
  if (n < k) return 0;
  T res = 1;
  REP(T, i, k) res = res * (n-i) % mod * modpow(i+1, mod-2) % mod;
  return res;
}

/*===========================================================================*/
//階乗とコンビネーションの計算(遅いバージョン)
VLL fac(10000);

template<class T>
T modfactorial(T n, T mod)
{
  /*
  n! を mod で割った余りを求める。
  グローバルに計算結果を格納する配列 VLL fac を用意する。初期値は -1。
  fac[n] は n! を mod で割った余り示す。
  */
  if (fac[n] != -1) return fac[n];
  REP(T, i, n+1)
  {
    if (fac[i] != -1) continue;
    if (i == 0) fac[i] = 1;
    else fac[i] = i * fac[i-1] % mod;
  }
  return fac[n];
}

template<class T>
T modfactorial2(T n, T mod)
{
  /*
  再帰関数で n! を mod で割った余りを求める。
  再帰回数が多くなると使えない。
  グローバルに計算結果を格納する配列 VLL fac を用意する。初期値は -1。
  fac[n] は n! を示す。
  */
  if (n == 0) return fac[0] = 1;
  if (fac[n] != -1) return fac[n];
  else return fac[n] = n * modfactorial2(n-1, mod) % mod;
}

template<class T>
T ModnCk(T n, T k, T mod)
{
  // フェルマーの小定理から nCk を mod で割った余りを求める。

  T res = modfactorial(n, mod);
  res = res * modpow(modfactorial(k, mod), mod-2, mod) % mod;
  res = res * modpow(modfactorial(n-k, mod), mod-2, mod) % mod;
  return res;
}

template<class T>
T modfact(T n, T mod)
{
  //n! を mod で割った余りを求める。
  T res = 1;
  REP2(T, i, 2, n+1)
  {
    res *= i;
    res %= mod;
  }
  return res;
}

/*============================================================================*/
//先に階乗のテーブルとその逆元テーブルを高速に構築し、 nCk (mod) を O(1) で求める。

const LL limit = 200000;  // n の上限値
VLL modfac(limit);        // n! を mod で割った余り
VLL invmodfac(limit);     // フェルマーの小定理における n! の逆元

template<class T>
void setmodfactorial(T _limit, T mod)
{
  REP(T, i, _limit+1)
  {
    if (i == 0) modfac[i] = 1;
    else modfac[i] = i * modfac[i-1] % mod;
  }
}

template<class T>
void setinvmodfactorial(T _limit, T mod)
{
  invmodfac[_limit] = modpow(modfac[_limit], mod-2, mod);
  REPR(T, i, _limit) invmodfac[i] = invmodfac[i+1] * (i+1) % mod;
}

template<class T>
T modnCk(T n, T k, T mod)
{
  return modfac[n] * invmodfac[k] % mod * invmodfac[n-k] % mod;
}

/*============================================================================*/
/*
コンストラクタで階乗とその逆元テーブルを構築。
テーブルを構築せずに modnCk2 を使いたければオブジェクトの宣言時に引数として mod のみを渡せばよい。
このクラスを宣言することで mod を法とした n! とその逆元、 a^b、 nCk の4つが求められる。
*/
template<class T>
struct ModCalc
{
  T _mod;
  V<T> modfac, invmodfac;

  ModCalc(const T& limit, const T& mod): _mod(mod), modfac(limit+1), invmodfac(limit+1)
  {
    setmodfactorial(limit+1);
    setinvmodfactorial(limit+1);
  }
  ModCalc(const T& mod): _mod(mod){}

  void setmodfactorial(const T& limit)
  {
    REP(T, i, limit)
    {
      if (i == 0 || i == 1) modfac[i] = 1;
      else modfac[i] = i * modfac[i-1] % _mod;
    }
  }

  void setinvmodfactorial(const T& limit)
  {
    invmodfac[limit-1] = modpow(modfac[limit-1], _mod-2);
    REPR(T, i, limit-1) invmodfac[i] = invmodfac[i+1] * (i+1) % _mod;
  }

  T modpow(T a, T b)
  {
    T res = 1;
    while (b > 0)
    {
      if (b & 1) res = res * a % _mod;
      a = a * a % _mod;
      b >>= 1;
    }
    return res;
  }

  // 1 <= k <= n <= 1e7
  T modnCk(const T& n, const T& k)
  {
    return modfac[n] * invmodfac[k] % _mod * invmodfac[n-k] % _mod;
  }

  // 1 <= n <= 1e9, 1 <= k <= 1e7
  T modnCk2(const T& n, const T& k)
  {
    if (n < k) return 0;
    T res = 1, tmp = 1;
    REP(T, i, k) res = res * (n-i) % _mod;
    REP(T, i, k) tmp = tmp * (i+1) % _mod;
    tmp = modpow(tmp, _mod-2);
    res = res * tmp % _mod;
    return res;
  }
};

/*============================================================================*/

template<class T>
void toAnotherMod(T &a, T &b, T &c, V<T>& rem_array)
{
  /*
  aで割ったときの余りがbである数字をcで割った余り(x = a*n + b とおいたとき、xをcで割った余り)を求める。
  結果を格納するための空配列を引数の最後に渡す。
  */
  int e_1 = a % c, e_2 = b % c;
  REP(T, i, c) rem_array.emplace_back((e_1 * i + e_2) % c);
  UNIQUE(rem_array);
}

/*============================================================================*/

// modをとる整数型
const int mod = 1e9+7;
struct mint
{
  LL x;

  mint (LL _x=0): x((_x%mod+mod)%mod){}

  mint operator-() const
  {
    mint tmp(-x);
    return tmp;
  }
  mint& operator+=(const mint& a)
  {
    if ((x += a.x) >= mod) x -= mod;
    return *this;
  }
  mint& operator-=(const mint& a)
  {
    if ((x += mod-a.x) >= mod) x -= mod;
    return *this;
  }
  mint& operator*=(const mint& a)
  {
    (x *= a.x) %= mod;
    return *this;
  }
  mint operator+(const mint& a) const
  {
    mint tmp(x);
    tmp += a.x;
    return tmp;
  }
  mint operator-(const mint& a) const
  {
    mint tmp(x);
    tmp -= a.x;
    return tmp;
  }
  mint operator*(const mint& a) const
  {
    mint tmp(x);
    tmp *= a.x;
    return tmp;
  }

  mint pow(LL b) const
  {
    mint res(1), a = *this;
    while (b > 0)
    {
      if (b & 1) res *= a;
      a *= a;
      b >>= 1;
    }
    return res;
  }

  mint inv() const{ return pow(mod-2); }
  mint& operator/=(const mint& a){ return *this *= a.inv(); }
  mint operator/(const mint& a) const
  {
    mint tmp(x);
    tmp /= a;
    return tmp;
  }

  bool operator>(const mint& a) const{ return (x > a.x); }
  bool operator<(const mint& a) const{ return (x < a.x); }
  bool operator==(const mint& a) const{ return (x == a.x); }
  bool operator>=(const mint& a) const{ return (x >= a.x); }
  bool operator<=(const mint& a) const{ return (x <= a.x); }
};
istream& operator>>(istream& is, mint& a){ return is >> a.x; }
ostream& operator<<(ostream& os, const mint& a){ return os << a.x; }

/*============================================================================*/

// mint想定
template<class T>
struct ModCalc
{
  V<mint> modfac, invmodfac;

  ModCalc(const T& limit): modfac(limit+1), invmodfac(limit+1)
  {
    setmodfactorial(limit+1);
    setinvmodfactorial(limit+1);
  }
  ModCalc(){}

  void setmodfactorial(const T& limit)
  {
    REP(T, i, limit)
    {
      if (i == 0 || i == 1) modfac[i] = 1;
      else modfac[i] = modfac[i-1] * i;
    }
  }

  void setinvmodfactorial(const T& limit)
  {
    invmodfac[limit-1] = mint(1) / modfac[limit-1];
    REPR(T, i, limit-1) invmodfac[i] = invmodfac[i+1] * (i+1);
  }

  // 1 <= k <= n <= 1e7
  mint modnCk(const T& n, const T& k)
  {
    return modfac[n] * invmodfac[k] * invmodfac[n-k];
  }

  // 1 <= n <= 1e9, 1 <= k <= 1e7
  mint modnCk2(const T& n, const T& k)
  {
    if (n < k) return 0;
    mint res = 1, tmp = 1;
    REP(T, i, k) res *= n-i;
    REP(T, i, k) tmp *= i+1;
    res /= tmp;
    return res;
  }
};

////////////////////////////////////////////////////////////////////////////////

// SPF(Smallest Prime Factor)を前計算して素因数分解する
template <class T>
struct PrimeFactor
{
  V<T> spf;

  PrimeFactor(const T& N): spf(N+1)
  {
    // SPFを求める。計算量はO(NloglogN)
    for (T i = 0; i <= N; i++) spf[i] = i;
    for (T i = 2; i*i <= N; i++) if (spf[i] == i) for(T j = i*i; j <= N; j += i) if (spf[j] == j) spf[j] = i;
  }

  // 素因数分解する。計算量はO(logn)
  map<T, T> prime_factorize(T n)
  {
    map<T, T> m;
    while(n > 1)
    {
      m[spf[n]]++;
      n /= spf[n];
    }
    return m;
  }

  // nの約数の個数を求める。計算量はO(logn)
  T get_divisor_num(const T& n)
  {
    map<T, T> pf = prime_factorize(n);
    T res = 1;
    for (auto p: pf) res *= p.second+1;
    return res;
  }
};

template<class T>
map<T, T> prime_factorize(T n)
{
  // 前処理なしで素因数分解する。計算量はO(√N)。
  map<T, T> res;
  for (T s = 2; n >= s * s; s += 2)
  {
    if (n % s == 0)
    {
      while(n % s == 0)
      {
        n /= s;
        res[s]++;
      }
    }
    if (s == 2) s--;
  }
  if (n != 1) res[n]++;
  return res;
}

template<class T>
map<T, T> pff(T N)
{
  //上の関数を用いて、N!を素因数分解した結果を返す。
  map<T, T> res;
  REP2(T, i, 2, N+1)
  {
    map<T, T> pre = prime_factorize(i);
    REPx(x, pre) res[x.first] += x.second;
  }
  return res;
}

template<class T>
bool primarity_test(T a)
{
  /*
  素数判定。
  a が素数なら true を返す。
  */
  if (a < 2) return false;
  else if (a == 2) return true;
  else if (a % 2 == 0) return false;
  for (T i = 3; i*i <= a; i += 2) if (a % i == 0) return false;
  return true;
}

template<class T>
struct Eratosthenes
{
  /*
  エラトステネスの篩で 0~limit までの素数を調べる。
  is_prime[i] が true なら i が素数であることを示す。
  計算量はO(nloglogn)。
  */
  VB is_prime;

  Eratosthenes(const T& limit): is_prime(limit+1, true)
  {
    for (T i = 0; i*i <= limit; ++i)
    {
      if (i == 0 || i == 1) is_prime[i] = false;
      else if (is_prime[i]) for (T j = i*i; j <= limit; j+=i) is_prime[j] = false;
    }
  }
};

template<class T>
T XOR_n(T n)
{
  // 0~nまでのXORを計算する。
  T count_1, a;
  if (n % 2)
  {
    count_1 = (n+1)/2;
    if (count_1 % 2) a = 1;
    else a = 0;
    return a;
  }
  else
  {
    count_1 = n/2;
    if (count_1 % 2) a = 1;
    else a = 0;
    return a ^ n;
  }
}

template<class T>
struct TopologicalSort
{
  VB visited, explored;
  V<T> Array; // ソート結果を逆順に格納する
  bool is_cyclic = false; // 閉路を検出した場合はtrueになる

  TopologicalSort(const VV<T>& G): visited(SIZE(G, T), false), explored(SIZE(G, T), false), Array(SIZE(G, T))
  {
    REP(T, i, SIZE(G, T))
    {
      if (sort(G, i))
      {
        is_cyclic = true;
        break;
      }
    }
  }

  bool sort(const VV<T>& G, const T& u)
  {
    static T tail = -1;
    if (visited[u]) return true;
    else if (explored[u]) return false;
    visited[u] = true;
    REPx(v, G[u]) if (sort(G, v)) return true;
    visited[u] = false;
    explored[u] = true;
    Array[++tail] = u;
    return false;
  }
};

template<class T>
struct BellmanFord
{
  V<T> dis;
  T flg;

  BellmanFord(const VV<P<T, T>>& G, const T& start): dis(SIZE(G, T), INF)
  {
    flg = bellman_ford(G, start);
  }

  bool bellman_ford(const VV<P<T, T>>& G, const T& s)
  {
    /*
    ベルマンフォード法。任意の始点から各点への最短距離を求める。計算量は O(VE)。
    引数のグラフは 0-indexed の隣接リストを想定。要素は、{to, weight}。
    始点から各点までの距離をdisに格納する。
    負の閉路を検出した場合はfalseを返し、正常に処理が終了した場合はtrueを返す。
    */
    dis[s] = 0;
    T V = SIZE(G, T);
    REP(T, roop, V) REP(T, i, V) REPx(x, G[i])
    {
      T u = i, v = x.first; // u -> v の辺
      T w = x.second;       // 辺の重み
      if (dis[v] > dis[u] + w)
      {
        dis[v] = dis[u] + w;
        if (roop == V-1) return false;
      }
    }
    return true;
  }
};

template<class T>
struct Dijkstra
{
  V<T> dis;

  Dijkstra(const VV<P<T, T>>& G, const T& start): dis(SIZE(G, T), INF)
  {
    dijkstra(G, start);
  }

  void dijkstra(const VV<P<T, T>>& G, const T& s)
  {
    /*
    ダイクストラ法。任意の始点から各点への最短距離を求める。計算量は O((V+E)logV)。
    引数のグラフは 0-indexed の隣接リストを想定。要素は、{to, weight}。
    始点から各点までの距離をdisに格納する。
    */
    dis[s] = 0;
    PRIORITY_QUEUE<P<T, T>> q; // 値の小さい要素を優先
    q.push({0, s}); // キューの要素は {weight, to}
    while(!q.empty())
    {
      T u = q.top().second;   // 頂点
      T c = q.top().first;    // 始点からのコスト
      q.pop();
      if (c > dis[u]) continue; // 枝刈り
      REPx(x, G[u])
      {
        T v = x.first;  // 行先
        T w = x.second; // 辺の重み
        if (c + w < dis[v])
        {
          dis[v] = c + w;
          q.push({dis[v], v});
        }
      }
    }
  }
};
/* ベルマンフォード、ダイクストラ法のサンプルコード
int main(){
  cin.tie(0);
  ios::sync_with_stdio(false);
  int N, M;
  cin >> N >> M;
  vector<VPII> G(N, VPII());
  REP(int, i, M){
    int u, v, l;
    cin >> u >> v >> l;
    u--; v--; // 1-indexed から 0-indexed へ
    G[u].emplace_back(v, l);
    G[v].emplace_back(u, l);
  }
  Dijkstra<int> D(G, 0); // BellmanFord
  debug1(D.dis);
  return 0;
}
input
------------
6 9
1 2 5
1 3 4
1 4 2
2 3 2
2 5 6
3 4 3
3 6 2
4 6 6
5 6 4
============
output
------------
0 5 4 2 10 6
*/

template<class T>
void warshall_floyd(VV<T>& G)
{
  /*
  ワ―シャルフロイド法。全点対間の最短距離を求める。計算量は O(V^3)。
  負の閉路がなければ負の辺にも対応可。
  引数のグラフは 0-indexed の隣接行列で与える必要がある。
  G[i][j] : 点iから点jへの辺の重み。辺が存在しない所はINFで初期化しておく。
  最終的にG[i][j]はiからjへの最小コストを格納する。
  中継点をkとしたとき、
  G[i][j] = G[i][k] + G[k][j]
  と表せることから、DPによりすべての点について、その点を中継点としたときのiからjへのパスの組み合わせをすべて計算することで、最終的な最小コストが求まる。
  */
  T V = SIZE(G, T);
  REP(T, i, V) G[i][i] = 0;
  REP(T, k, V) REP(T, i, V) REP(T, j, V) G[i][j] = min(G[i][j], G[i][k] + G[k][j]);
}
/*ワーシャルフロイド法のサンプルコード
int main(){
  cin.tie(0);
  ios::sync_with_stdio(false);
  int N, M;
  cin >> N >> M;
  VVI G(N, VI(N, INF));
  REP(int, i, M)  {
    int u, v, l;
    cin >> u >> v >> l;
    u--; v--;                   // 1-indexed から 0-indexed へ
    G[u][v] = l;
    G[v][u] = l;
  }
  warshall_floyd(G);
  ENDL
  REP(int, i, N) REP(j, N) cout << i << "から" << j << "へのコスト : " << G[i][j] << endl;
  return 0;
}
*/

int h, w;
template<class T1, class T2>
T1 bfs(VV<T2>& MAP, T1 sy, T1 sx, T1 gy, T1 gx)
{
  // 迷路探索のテンプレ
  const static V<T1> px = {0, 0, 1, -1}, py = {1, -1, 0, 0};
  VV<T1> dis(h, V<T1>(w, -1));
  queue<P<T1, T1>> q;
  q.push({sx, sy});
  dis[sy][sx] = 0;
  while(!q.empty())
  {
    T1 y = q.front().second, x = q.front().first;
    q.pop();
    REP(T1, i, 4)
    {
      T1 yy = y + py[i], xx = x + px[i];
      if (0 <= xx && xx < w && 0 <= yy && yy < h && MAP[yy][xx] != '#' && dis[yy][xx] == -1)
      {
        if (yy == gy && xx == gx) return dis[y][x] + 1;
        dis[yy][xx] = dis[y][x] + 1;
        q.push({xx, yy});
      }
    }
  }
  return -1;
}

template<class T1, class T2>
T1 bfs01(VV<T2>& M, T1 sy, T1 sx, T1 gy, T1 gx)
{
  /*
  01-BFSのテンプレ(迷路)
  dequeを使うことで、辺長が 0 or 1 のグラフ最短経路問題を解く
  BFSもダイクストラ法も、始点からの距離が暫定的に最短の頂点を順に取り出して、その点から伸びる辺でその先の頂点の暫定最短距離を更新することを繰り返す
  このアルゴリズムは、保持しているデータから常に最小のものを取り出せるようなデータ構造があれば実現できる
  迷路探索のような通常のBFSの場合は、辺長が1の単一始点グラフ最短経路問題を解くことなので始点から順に探索保持していき、データを取り出すときは入れた順に取り出せばよいのでqueueを使えばよい(FIFO)
  ダイクストラ法の場合、単一始点非負辺長グラフ最短経路問題を解くことなので、まず始点から探索を始め、始点からの暫定最短距離が最小の頂点vを取り出し、vから伸びる辺によって暫定最短距離を更新できる頂点があれば、更新してその点を保持していけばよい。これはpriority_queueによって実現できる
  01-BFSとは、辺長が0または1の単一始点グラフ最短経路問題を解くようなアルゴリズム。通常のBFSと同様に、始点から探索をはじめ、今見ている頂点から伸びた先の頂点が未更新の場合は順に更新して保持していけばよい。しかしqueueによって保持していくとデータを取り出すときに辺長が1だった頂点と辺長が0だった点が逆順に取り出され得る。そこでpriority_queueを使ったダイクストラ法で解くこともできるが、dequeによってさらに効率的に解くことができる。
  具体的には、頂点の暫定最短距離の更新の際、dequeへデータを入れるには、0の辺を使う場合は先頭に、1の辺を使う場合は末尾に加えればよい。そうすればコストの大きい頂点が後ろ、小さい頂点が前の方へ格納されることになり、常にdequeの先頭からデータを取り出すようにすることで所望のデータ構造が得られる。
  */
  static V<T1> px = {0, 0, 1, -1}, py = {1, -1, 0, 0};
  VV<T1> dis(h, V<T1>(w, -1));
  dis[sy][sx] = 0;
  deque<P<T1, T1>> q;
  q.push_back({sx, sy});
  while(!q.empty())
  {
    T1 x = q.front().first, y = q.front().second;
    q.pop_front();
    REP(T1, i, 4)
    {
      T1 xx = x + px[i], yy = y + py[i];
      if(0 <= xx && xx < w && 0 <= yy && yy < h && dis[yy][xx] == -1)
      {
        if (M[yy][xx] == '#')
        {
          dis[yy][xx] = dis[y][x] + 1;
          q.push_back({xx, yy});
        }
        else
        {
          dis[yy][xx] = dis[y][x];
          q.push_front({xx, yy});
        }
        if (yy == gy && xx == gx) return dis[yy][xx];
      }
    }
  }
}

template<class T>
bool inline is_colorable(VV<T>& G, T now, T n, T v)
{
  /*
  グラフの頂点に 0 ~ n-1 までの n 種類の番号を付けるとき、隣り合う頂点同士の番号の差がちょうど 1 もしくは n-1 になるようにできるかを判定する。(n-1 の次の番号は 0 とする) n = 2 のときは二部グラフの判定をすることになる。
  グラフは 0-indexed の隣接リストを想定。 now には 0 を引数として渡し、0番目の頂点から探索を開始するようにする。 v には頂点の数を渡す。
  */
  static V<T> color(v, -1);
  if (now == 0) color[now] = 0;
  REPx(x, G[now])
  {
    if (color[x] != -1)
    {
      if (abs(color[now]-color[x]) != 1 && abs(color[now]-color[x]) != n-1) return false;
      continue;
    }
    color[x] = (color[now] + 1) % n;
    if (!is_colorable(x, n, v)) return false;
  }
  return true;
}

template<class T>
bool inline is_colorable_a(VV<T>& G, T now, T n, T v)
{
  //上の隣接行列バージョン
  static V<T> color(v, -1);
  if (now == 0) color[now] = 0;
  REP(T, next, v)
  {
    if (G[now][next] == 1)
    {
      if (color[next] != -1)
      {
        if (abs(color[now]-color[next]) != 1 && abs(color[now]-color[next]) != n-1) return false;
        continue;
      }
      color[next] = (color[now] + 1) % n;
      if (!is_colorable_a(next, n, v)) return false;
    }
  }
  return true;
}

template<class T>
int LIS(V<T> A)
{
  // 数列 A のLIS(Longest Increasing Subsequence)の長さを求める。計算量はO(NlogN)
  int N = SIZE(A, int);
  V<T> dp(N, INF);
  REP(int, i, N) dp[bisect_left(dp, A[i])] = A[i];
  return bisect_left(dp, INF);
}

int LCS(auto A, auto B)
{
  // 配列 A と B のLCS(Longest Common Subsequence)の長さを求める。計算量はO(NM)
  int N = SIZE(A, int), M = SIZE(B, int);
  VVI dp(N+1, VI(M+1, 0));
  REP(int, i, N) REP(int, j, M)
  {
    dp[i+1][j+1] = max(dp[i+1][j], dp[i][j+1]);
    if (A[i] == B[j]) dp[i+1][j+1] = dp[i][j] + 1;
  }
  return dp[N][M];
}

string D2NC(ULL D, int N)
{
  // 10進数表記の数字 D をN進数表記に変換する
  if (N <= 1)
  {
    cout << "Nは2以上の整数を入力してください" << endl;
    exit(0);
  }
  string res;
  string alp = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  while (D > 0)
  {
    if (D % N > 9) res += alp[D % N - 10];
    else res += '0' + D % N;
    D /= N;
  }
  if (res.empty()) res = "0";
  REVERSE(res);
  return res;
}

ULL POW(ULL a, ULL b)
{
  // a^b を計算する
  ULL res = 1;
  REP(int, i, b) res *= a;
  return res;
}

ULL N2DC(string V, ULL N)
{
  // N進数表記の数字 V を10進数表記に変換する
  map<char, int> check;
  string number = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  REP(int, i, 36) check[number[i]] = i;
  REPx(&x, V)
  {
    if (islower(x)) x -= 32;
    if (check[x] >= (int)N)
    {
      cout << V << "は" << N << "進数ではありません" << endl;
      exit(0);
    }
  }
  ULL res = 0;
  REP(ULL, i, SIZE(V, ULL))
  {
    ULL E = SIZE(V, ULL) - i - 1;
    if (isdigit(V[i])) res += (V[i]-'0') * POW(N, E);
    else
    {
      char x = V[i];
      res += (x - 'A' + 10) * POW(N, E);
    }
  }
  return res;
}

string M2NC(string V, ULL M, ULL N)
{
  // M進数表記の数字 V をN進数表記に変換する
  // ULLを超えた数字は扱えないのでpythonで計算したほうが安全
  if (N <= 1)
  {
    cout << "Nは2以上の整数を入力してください" << endl;
    exit(0);
  }
  map<char, int> check;
  string number = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  REP(int, i, 36) check[number[i]] = i;
  REPx(&x, V)
  {
    if (islower(x)) x -= 32;
    if (check[x] >= (int)M)
    {
      cout << "入力された数字は" << M << "進数ではありません" << endl;
      exit(0);
    }
  }
  if (M == N) return V;
  if (M == 10) return D2NC(stoull(V), N);
  if (N == 10) return to_string(N2DC(V, M));
  return D2NC(N2DC(V, M), N);
}

/*
void calcPI(int &n)
{
  //ガウス・ルジャンドルのアルゴリズムによって円周率を求める。小数第m位程度の精度が欲しいときは、log2(m)回程度の反復計算を行えばよい。つまり、n = log2(m)を引数に入れる。
  //グローバルに以下の三行を追加する必要がある。※コンパイルに時間がかかる
  //    #include <boost/multiprecision/cpp_dec_float.hpp>
  //    namespace mp = boost::multiprecision;
  //    using cpp_dec_float_1000 = mp::number<mp::cpp_dec_float<1000> >;
  cpp_dec_float_1000 x = 2.0f;
  cpp_dec_float_1000 a = 1.0, b = 1 / mp::sqrt(x), t = 0.25, p = 1.0;
  cpp_dec_float_1000 an, bn, tn, pn;
  REP(i, n)
  {
    an = (a + b) / 2;
    bn = mp::sqrt(a * b);
    tn = t - p * (a - an) * (a - an);
    pn = p * 2;
    a = an; b = bn; t = tn; p = pn;
  }
  cpp_dec_float_1000 res = (a + b) * (a + b) / (4 * t);
  cout << setprecision(numeric_limits<decltype(a)>::digits10 + 1) << res << endl;
}
*/

// 一次元累積和を計算して、指定された範囲内の総和を求める。
template<class T>
struct Cumulative_Sum
{
  V<T> s;
  T N;

  Cumulative_Sum(const V<T>& a): N(SIZE(a, T)), s(N+1)
  {
    /*
    Parameters
    ----------
    a : 累積和を求めたい数列
    N : aのサイズ
    s : 累積和計算用テーブル
    */
    REP(int, i, N) s[i+1] = s[i] + a[i];
  }

  T calc_sum(const T& l, const T& r)
  {
    /*
    Parameters
    ----------
    l, r : 総和を求めたい範囲。1-indexedで[l, r]の範囲を示す。

    Returns
    -------
    [l, r]内の総和
    */
    return s[r] - s[l-1];
  }
};
/*
-------------
N
a1 a2 ... aN
Q
l1 r1
l2 r2
.  .
.  .
.  .
lQ rQ
-------------
1-indexedの数列を想定している。
l, r : [l, r]の範囲を表す。
Q個のクエリで指定された範囲内の総和を求める。
===========
<input>
5
1 3 2 5 4
3
2 4
1 5
1 3
===========
<output>
10
15
6
===========
int main(){
  int N, Q;
  cin >> N;
  VI A(N);
  INPUT(A);
  cin >> Q;
  Cumulative_Sum<int> cs(A);
  REP(int, i, Q)  {
    int l, r;
    cin >> l >> r;
    cout << cs.calc_sum(l, r) << endl;
  }
  return 0;
}
*/

// 二次元累積和を計算して、指定された長方形領域内の総和を求める。
template<class T>
struct Cumulative_Sum_2d
{
  T H, W;
  VV<T> s;
  Cumulative_Sum_2d(const VV<T>& a): H(SIZE(a, T)), W(SIZE(a[0], T)), s(H+1, V<T>(W+1))
  {
    /*
    Parameters
    ----------
    a : 累積和を求めたい二次元配列
    H : aの高さ
    W : aの幅
    s : 累積和計算用テーブル。s[y][x]は、1-indexedで[1, y]×[1, x]
             領域(原点を含む高さy、幅xの長方形領域)内の総和を表す。
    */
    REP(T, y, H) REP(T, x, W) s[y+1][x+1] = s[y][x+1] + s[y+1][x] - s[y][x] + a[y][x];
  }

  T calc_sum(const T& yh, const T& yt, const T& xh, const T& xt)
  {
    /*
    Parameters
    ----------
    yh, yt : 総和を求めたい長方形領域の縦の範囲。1-indexedで[yh, yt]
             の範囲を示す。
    xh, xt : 総和を求めたい長方形領域の横の範囲。1-indexedで[xh, xt]
             の範囲を示す。

    Returns
    -------
    [yh, yt]×[xh, xt]領域内の総和
    */
    return s[yt][xt] - s[yt][xh-1] - s[yh-1][xt] + s[yh-1][xh-1];
  }
};
/*
----------------
H W
a11 a12 ... a1W
 .   .  ...  .
 .   .  ...  .
 .   .  ...  .
aH1 aH2 ... aHW
Q
yh1 yt1 xh1 xt1
 .   .   .   .
 .   .   .   .
 .   .   .   .
yhQ ytQ xhQ xtQ
----------------
1-indexedの二次元領域を想定している。
yh yt xh xt : [yh, yt]×[xh, xt]の長方形領域を表す。
Q個のクエリで指定された領域内の総和を求める。
===========
<input>
4 5
1 8 7 3 2
9 1 3 4 6
3 5 8 1 4
2 7 3 2 5
3
2 3 3 5
1 2 2 3
1 4 1 5
===========
<output>
26
19
84
===========
int main(){
  int H, W, Q;
  cin >> H >> W;
  VVI A(H, VI(W));
  INPUT2(A);
  cin >> Q;
  Cumulative_Sum_2d<int> cs2(A);
  REP(int, i, Q)  {
    int yh, yt, xh, xt;
    cin >> yh >> yt >> xh >> xt;
    cout << cs2.calc_sum(yh, yt, xh, xt) << endl;
  }
  return 0;
}
*/

// フィボナッチ数列の第n項目までを計算する
template<class T>
struct Fibonacci
{
  V<T> memo;
  Fibonacci(T n): memo(n+1, 0){}

  T calc(T m)
  {
    if (m == 1 || m == 2) return 1;
    if (memo[m]) return memo[m];
    return memo[m] = calc(m-1) + calc(m-2);
  }
};

template<class T>
struct UnionFind
{
  V<T> parent;
  UnionFind(const T& n): parent(n, -1){}

  // xの親を検索
  T find(const T& x)
  {
    if (parent[x] < 0) return x;
    return parent[x] = find(parent[x]);
  }

  // 同じ集合に属するかどうかをチェック
  bool same_check(const T& x, const T& y)
  {
    return find(x) == find(y);
  }

  // 自分の属する集合のサイズを返す
  T size(const T& x)
  {
    T y = find(x);
    return -parent[y];
  }

  // xとyの属する集合を併合
  bool unite(const T& x, const T& y)
  {
    T x_root = find(x), y_root = find(y);
    if (x_root == y_root) return false;
    else
    {
      if (-parent[x_root] < -parent[y_root]) swap(x_root, y_root);
      parent[x_root] += parent[y_root];
      parent[y_root] = x_root;
      return true;
    }
  }
};

int main()
{
  cin.tie(0);
  ios::sync_with_stdio(false);
  int n;
  cin >> n;
  cout << XOR_n(n);
}