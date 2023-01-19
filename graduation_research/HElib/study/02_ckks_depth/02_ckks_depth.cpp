#include <helib/helib.h>

using namespace std;
using namespace helib;

/* 
以下のメモは下記のソースコードを要約したものである。
https://github.com/homenc/HElib/blob/master/examples/tutorial/02_ckks_depth.cpp

概要
+ CKKS方式では準同型演算を繰り返すことで、暗号文にノイズ("noise")が溜まっていく。
+ 準同型演算の計算回数を深さ("depth")と呼び、これが増加するほどノイズが増え、容量("capacity")と精度("accuracy")が低下する。

容量について
+ パラメータ"bits"より少し小さい値から始まる
+ 準同型演算をするごとに減っていき、1を下回ると復号できなくなる

精度
+ パラメータ"precision"の値をもとに、約 2^{-precision} より大きくならないようにする必要がある。
+ 暗号文と平文の間に生まれる誤差の程度を示している
+ 誤差は max|c[i] - p[i]| (i = 1, ..., n - 1) で表される(n = slots)
+ 準同型演算をするごとに誤差が増加していく

2種のパラメータの確認方法
容量：chipertext.capacity()
精度：chipertext.errorBound() (この場合は先の式で算出される誤差を示す)
*/

int main() {
    Context context = ContextBuilder<CKKS>()
    .m(32 * 1024)
    .bits(358)//358
    .precision(16)
    .c(6).build();

    cout << "securityLevel=" << context.securityLevel() << "\n";
    
    long n = context.getNSlots();
    
    SecKey secretKey(context);
    secretKey.GenSecKey();
    const PubKey& publicKey = secretKey;

    vector<double> v(n);
    for (long i = 0; i < n; ++i) {
        v[i] = sin(2.0 * PI * i / n);
    }

    PtxtArray p(context, v);
    Ctxt c(publicKey);
    p.encrypt(c);


    cout << "c.capacity=" << c.capacity() << " ";
    cout << "c.errorBound=" << c.errorBound() << "\n";

    for(int i = 0; i < 5; ++i) {
        c *= c;
        p *= p;

        cout << "c.capacity=" << c.capacity() << " ";
        cout << "c.errorBound=" << c.errorBound() << "\n";
    }

    PtxtArray pp(context);
    pp.decrypt(c, secretKey);

     double distance = Distance(p, pp);
    cout << "distance=" << distance << "\n";
}