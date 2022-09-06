#include <iostream>
#include <helib/helib.h>
#include <chrono>
using namespace std;
using ll = long long;

// (p = 7, m = 12351), (p = 17, m = 14351)
const unsigned long p = 7;
unsigned long m = 12351;
unsigned long r = 1;
unsigned long bits = 500;
unsigned long c = 2;

vector<vector<ll>> ConstDivCoefficients(p, vector<ll>(p));
vector<vector<ll>> ConstEqCoefficients(p, vector<ll>(p));

// Pows関数ではべき乗列の作成を工夫することで、乗算が深くなることを防いでいる
vector<helib::Ctxt> Pows(helib::Ctxt c) {
	int l = log2(p);
	vector<helib::Ctxt> c_pows;

	// c_powsをp個分の配列にresizeする
	for (int i = 0; i < p - 1; ++i) c_pows.push_back(c);

	// c^1からc^l-1までを計算する
	for (int i = 0; i < l; ++i) {
		int base = pow(2, i);

		for (int j = 0; j < base; ++j) {
			helib::Ctxt temp(c_pows[base - 1]);
			temp.multiplyBy(c_pows[j]);
			c_pows[base + j] = temp;
		}
	}

	// c^lからc^p-1までを計算する
	if (pow(2, l) < p - 1) {
		long long base = pow(2, l);

		for (int i = 0; i < p - base - 1; ++i) {
			helib::Ctxt temp(c_pows[base - 1]);
			temp.multiplyBy(c_pows[i]);
			c_pows[base + i] = temp;
		}
	}

	return c_pows;
}

// lagrange_interpolation関数で使われる
ll modpow(ll x, ll n, const ll mod) {
	ll res = 1;

	while (n > 0) {
		if (n & 1) (res *= x) %= mod;
		(x *= x) %= mod;
		n >>= 1;
	}
	return res;
}

// lagrange_interpolation関数で使われる
ll inverse(ll x, const ll mod) {
	return modpow(x, mod - 2, mod);
}

// ラグランジュ補間を行い係数を求める
// https://suikaba.hatenablog.com/entry/2019/08/11/021048
vector<ll> lagrange_interpolation(vector<ll> x, vector<ll> y) {
	const int n = x.size();

	for (int i = 0; i < n; ++i) {
		x[i] %= p;
		y[i] %= p;
	}

	// 係数c_0, c_1, ..., c_{n+1}の配列
	vector<ll> all_coefficient(n + 1);
	all_coefficient[0] = 1;

	// Π_{k = 0}^n (x - x_k)を計算し、係数c_0, c_1, ..., c_{n+1}を求める
	// (x - x0), (x - x1), ..., (x - xn)と順々に掛けていき、係数を計算する 
	for (int i = 0; i < n; ++i) {
		vector<ll> nxt(n + 1);

		// 新たに(x - x_i)を掛ける準備として、次数を1つ上げておく
		for (int k = 0; k < n; ++k) {
			nxt[k + 1] = all_coefficient[k];
		}

		// (x - x_i)を掛ける
		for (int k = 0; k < n; ++k) {
			nxt[k] = (p + nxt[k] - x[i] * all_coefficient[k] % p) % p;
		}

		all_coefficient = move(nxt);
	}

	vector<ll> coefficient(n);

	// Q_i = y_i / Π_{k != k}(x_i - x_k) を計算する
	// i: データ点(x0, y0), (x1, y1), ..., (xi, yi)
	// 各データ点について係数を計算していく
	for (int i = 0; i < n; ++i) {
		ll qi = 1;

		// Q_iの分母部分を計算する
		for (int k = 0; k < n; ++k) {
			if (i == k) continue; // Π_{k != i}のため、i == kの場合は計算しない
			qi = qi * (p + x[i] - x[k]) % p; // Π_{i != k}(x_i - x_k)を計算する
		}
		qi = y[i] * inverse(qi, p) % p; // y[i]とΠ_{i != k}(x_i - x_k)の逆数を掛けることでQ_iを求める


		// 多項式の係数c0, c1, ..., cnを求める
		vector<ll> temp_coefficient = all_coefficient;
		for (int k = n - 1; k >= 0; --k) {
			// Q_i と  Π_{k = i}^n (x - x_k) / (x - xi) を掛けて係数を算出する
			// c_k = q_i * Π_{i != k}(x - x_k)
			coefficient[k] = (coefficient[k] + qi * temp_coefficient[k + 1]) % p;

			// Π_{"k = i"}^n (x - x_k)上で考えるため、(x - x_i)で割る必要がある
			// その操作をこの計算で行う
			// c_j' = c_{k + 1} + c_{k+1}' * x_i
			temp_coefficient[k] = (temp_coefficient[k] + temp_coefficient[k + 1] * x[i]) % p;
		}
	}
	return coefficient;
}


// g(x) := x/d とし、f(x) = g(x)を満たす関数f(x)について多項式補完を行い、その係数を求める
// d in Z_p について計算する。
void PrecomupteConstDivCoefficient() {
	for(int d = 1; d < p; ++d) {
		vector<ll> data_x(p);
		vector<ll> data_y(p);

		for(int x = 0; x < p; ++x) {
			data_x[x] = x;
			data_y[x] = x / d;
		}

		ConstDivCoefficients[d] = lagrange_interpolation(data_x, data_y);
	}
}

// g(x) := 1 :if(x=y), 0 :if(x!=y) とし、f(x) = g(x)を満たす関数f(x)について多項式補完を行い、その係数を求める
// d in Z_p について計算する。
void PrecomupteConstEqCoefficient() {
	for(int d = 1; d < p; ++d) {
		vector<ll> data_x(p);
		vector<ll> data_y(p);

		for(int x = 0; x < p; ++x) {
			data_x[x] = x;
            if(x == d) data_y[x] = 1;
            else data_y[x] = 0;
		}

		ConstEqCoefficients[d] = lagrange_interpolation(data_x, data_y);
	}
}

helib::Ctxt ConstDiv(vector<helib::Ctxt> C_pow, int divisor, const helib::Context& context, const helib::PubKey& pk) {
	// C_div_result = C_0 を作成する
	helib::Ptxt<helib::BGV> P_div_result(context, vector<int>(1, 0));
	helib::Ctxt C_div_result(pk);
	pk.Encrypt(C_div_result, P_div_result);

	for(int i = 0; i < p; ++i) {
		if(i == 0) {
			C_div_result.addConstant(NTL::ZZX(ConstDivCoefficients[divisor][i]));
			continue;
		}

		helib::Ctxt temp(C_pow[i - 1]); // C_powはC_pow[0] = C^1, C_pow[1] = C^2と、添え字とべき数がズレているため、i - 1 している
		temp.multByConstant(NTL::ZZX(ConstDivCoefficients[divisor][i]));
		C_div_result.addCtxt(temp);
	}

	return C_div_result;
}

helib::Ctxt ConstEq(vector<helib::Ctxt> C_pow, int divisor, const helib::Context& context, const helib::PubKey& pk) {
	// C_div_result = C_0 を作成する
	helib::Ptxt<helib::BGV> P_div_result(context, vector<int>(1, 0));
	helib::Ctxt C_div_result(pk);
	pk.Encrypt(C_div_result, P_div_result);

	for(int i = 0; i < p; ++i) {
		if(i == 0) {
			C_div_result.addConstant(NTL::ZZX(ConstEqCoefficients[divisor][i]));
			continue;
		}

		helib::Ctxt temp(C_pow[i - 1]); // C_powはC_pow[0] = C^1, C_pow[1] = C^2と、添え字とべき数がズレているため、i - 1 している
		temp.multByConstant(NTL::ZZX(ConstEqCoefficients[divisor][i]));
		C_div_result.addCtxt(temp);
	}

	return C_div_result;
}

helib::Ctxt Div(helib::Ctxt C_a, helib::Ctxt C_d, const helib::Context& context, const helib::PubKey& pk) {
	// C_sum = C_0 を作成する
	helib::Ptxt<helib::BGV> P_sum(context, vector<int>(1, 0));
	helib::Ctxt C_sum(pk);
	pk.Encrypt(C_sum, P_sum);

	// 各暗号文のべき乗列を計算する
	vector<helib::Ctxt> C_a_pow = Pows(C_a);
	vector<helib::Ctxt> C_d_pow = Pows(C_d);

	for(int i = 1; i < p; ++i) {
		// C_a/i の準同型演算の結果を返す
		helib::Ctxt C_div_a_i = ConstDiv(C_a_pow, i, context, pk);

		// d == i の場合：C_0
		// d != i の場合：C_1
		helib::Ctxt C_equal_d_i = ConstEq(C_d_pow, i, context, pk);

		// C_sum = C_sum + FHE.Mult(C_a/i, C_d=i)
		// d == i の場合：C_sum = C_sum + C_a/i * C_0 = C_sum + C_0  = C_sum (C_0より0を足しているだけなので、C_sumのまま)
		// d != i の場合：C_sum = C_sum + C_a/i * C_1 = C_sum + C_a/i        (C_a/iの結果がC_sumに加算される)
		C_div_a_i.multiplyBy(C_equal_d_i);
		C_sum.addCtxt(C_div_a_i);
	}

	helib::Ctxt C_div_a_d(C_sum);

    return C_div_a_d;
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		cout << "few argument error" << endl;
		return 1;
	}

    // セットアップ
	helib::Context context = helib::ContextBuilder<helib::BGV>()
		.m(m)
		.p(p)
		.r(r)
		.bits(bits)
		.c(c)
		.build();

    // 秘密鍵
	helib::SecKey secret_key(context);
	secret_key.GenSecKey();
	helib::addSome1DMatrices(secret_key);

    // 公開鍵
	const helib::PubKey& public_key = secret_key;
	const helib::EncryptedArray& ea = context.getEA();

    // 平文
    int a = *argv[1] - '0', d = *argv[2] - '0';
	helib::Ptxt<helib::BGV> plaintext_a(context, vector<int>(1, a));
	helib::Ptxt<helib::BGV> plaintext_d(context, vector<int>(1, d));

    // 暗号化
	helib::Ctxt cihpertext_a(public_key);
	helib::Ctxt cihpertext_d(public_key);
	public_key.Encrypt(cihpertext_a, plaintext_a);
	public_key.Encrypt(cihpertext_d, plaintext_d);

	PrecomupteConstDivCoefficient();
	PrecomupteConstEqCoefficient();

	// 時間計測
	chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    helib::Ctxt cihpertext_div_a_d = Div(cihpertext_a, cihpertext_d, context, public_key);

	end = chrono::system_clock::now();

	double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);

	// 復号
	helib::Ptxt<helib::BGV> result(context);
	secret_key.Decrypt(result, cihpertext_div_a_d);
	cout << result << endl;
	printf("time %lf[ms]\n", time);
}
