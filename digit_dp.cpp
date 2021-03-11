#include <iostream>
#include <vector>

using namespace std;

template<class T> using V = vector<T>;
template<class T> using VV = V<V<T>>;
template<class T> using VVV = V<VV<T>>;
using LL = long long;
using VLL = V<LL>;
using VVLL = V<VLL>;

#define REP(type, i, n) for (type i = 0; i < (type)(n); ++i)
#define SIZE(a, type) ((type)(a).size())

const LL MOD = 1e9 + 7;

// N以下の整数を数える(0~NまでのN+1通り)
int main()
{
  cin.tie(0);
  ios::sync_with_stdio(false);
  string N;
  cin >> N;
  int digit = SIZE(N, int);  // Nの桁数
  VVLL dp(digit, VLL(2, 0));
  dp[0][1] = N[0] - '0';  // 0~N[0]-1 までのN[0]個は未満フラグが1
  dp[0][0] = 1;           // 0桁目がN[0]なら未満フラグは0
  REP(int, i, digit-1){
    const int D = N[i+1] - '0'; // i+1桁目の数字
    REP(int, j, 2){
      REP(int, d, (j ? 10 : D+1)){        // 未満フラグが立っているときは0~9まですべてが可能で、そうでないときは0~Dまでのみ可能
        dp[i+1][j || d < D] += dp[i][j];  // 遷移先の未満フラグは遷移前と同じで、d < Dのときは未満フラグを立てる
      }
    }
  }
  cout << dp[digit-1][0] + dp[digit-1][1] << endl;

  return 0;
}
/*
string N として、 N[i] はNのi桁目の数字を表す文字
dp[i][j] : 上位からi桁目までの数字の総数
           jは構成される数字がN未満であるかどうかを示すフラグ(未満フラグ)
            ・ j = 0 のとき構成される数字がi桁目まで完全に一致していることを示す
            ・ j = 1 のとき構成される数字がN未満となることを示す(0≦k≦iについてNのk桁目より小さくなるものが存在する)

(遷移)
int D = N[i+1]-'0';
dp[i][j]が確定しているものとして状態の遷移を考える
 ・ dp[i][1] のとき、すでにN未満の数字であることが確定しているので、i+1桁目の数字として考えられるのは、0~9までの数字すべてである
 ・ dp[i][0] のとき、i桁目の数字まですべて一致しているので、i+1桁目の数字として考えられるのは、0~Dまでの数字のみである(D以上の数字だとNを超えてしまう)
これらをまとめて考えると、その遷移の仕方は、d = (i+1桁目として選ぶ数字) として
 ・ dp[i][1] -> dp[i+1][1]
 ・ dp[i][0] -> dp[i+1][0] (if (d==D)) or dp[i+1][1] (if (d < D))
となる
これらを場合分けで実装してもいいが、よく見ると遷移前の未満フラグが立っているときはときは遷移先のフラグも立っていて、立っていないときは d < D のときに遷移先のフラグが立つことから実装を工夫している

(初期条件)
dp[0][0] = 1 (0桁目の数字まで一致している)
dp[0][1] = N[0]-'0' (Nの0桁目より小さい数字の数)

(最終状態)
int digit = SIZE(N);
dp[digit-1][0] : digit-1 桁目まででN未満の数字の個数(0以上N未満の整数の個数)
dp[digit-1][1] : digit-1 桁目まですべての数字がNと一致している数字の個数(Nの1つのみ)

(参照)
https://torus711.hatenablog.com/entry/20150423/1429794075
https://qiita.com/murai_mart/items/1a7a4d10abc0c3b5b53f
*/

//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////

// N以下の整数でいずれかの桁に3を含むものの総数を数える
int Main()
{
  cin.tie(0);
  ios::sync_with_stdio(false);
  string N;
  cin >> N;
  int digit = SIZE(N, int);
  VVV<LL> dp(digit, VVLL(2, VLL(2, 0)));
  REP(int, d, N[0]-'0'+1) dp[0][d < N[0]-'0'][d == 3]++;
  REP(int, i, digit-1)
  {
    int D = N[i+1] - '0';
    REP(int, j, 2) REP(int, k, 2) REP(int, d, j ? 10 : D+1) dp[i+1][j || d < D][k || d == 3] += dp[i][j][k];
  }
  cout << dp[digit-1][0][1] + dp[digit-1][1][1] << endl;

  return 0;
}
/*
dp[i][j][k] : 上位からi桁目までの数字で3を含むものの総数
              jは未満フラグ
              kは構成される数字に3が含まれているかどうかを示すフラグ
               ・ k = 0 のとき構成される数字に3が含まれないことを示す
               ・ k = 1 のとき構成される数字に3が含まれることを示す(0≦k≦iについてNのk桁目より小さくなるものが存在する)

(遷移)
int D = N[i+1] - '0';
・ dp[i][0][0] -> dp[i+1][0][0] (d == D != 3)
                 dp[i+1][0][1] (d == D == 3)
                 dp[i+1][1][0] (d != 3)
                 dp[i+1][1][1] (d == 3)
・ dp[i][0][1] -> dp[i+1][0][1] (d == D)
                 dp[i+1][1][1] (d < D)
・ dp[i][1][0] -> dp[i+1][1][0] (d != 3)
                 dp[i+1][1][1] (d == 3)
・ dp[i][1][1] -> dp[i+1][1][1] (for all d)
*/

//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////

// N以下の正整数であって、十進法表記したときの各桁の数の和がDの倍数であるものの個数を mod 1e9+7 で求める
LL f(string X, int s)
{
  int digit = SIZE(X, int);
  VVV<LL> dp(digit+1, VVLL(2, VLL(s, 0)));
  dp[0][0][0] = 1;
  REP(int, i, digit)
  {
    int D = X[i] - '0';
    REP(int, j, 2) REP(int, k, s) REP(int, d, j ? 10 : D+1)
    {
      dp[i+1][j || d < D][(k + d) % s] += dp[i][j][k];
      dp[i+1][j || d < D][(k + d) % s] %= MOD;
    }
  }
  return (dp[digit][0][0] + dp[digit][1][0] - 1) % MOD;
}
/*
dp[i][j][k] : 上位からi桁目までの数字で各桁の和の mod s がkになるものの総数
              jは未満フラグ
これまでと違って0桁はどの数字も選んでいない状態を示す(これまでと同様のやり方で実装してもいいが0桁から始めることで、dp[0][0][0] = 1 とし、以降の桁についてはdpテーブルを回すだけでよくなる)

(遷移)
int D = N[i+1] - '0';
・ dp[i][0][k] -> dp[i+1][0][(k+d)%s] (d == D)
                 dp[i+1][1][(k+d)%s] (d < D)
・ dp[i][1][k] -> dp[i+1][1][(k+d)%s] (for all d)
*/