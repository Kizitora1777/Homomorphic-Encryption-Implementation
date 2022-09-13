#include "openfhe.h"

using namespace std;
using namespace lbcrypto;

using ll = long long;

CryptoContext<DCRTPoly> cryptoContext;
KeyPair<DCRTPoly> keyPair;
const int p = 17;
// Okada, 2019, p.9 より、Pows関数の深さはlog2(p-1)となる
// 加えて、Div関数でも準同型乗算を行うので、深さを1つ加える
const int multiplicative_depth = log2(p - 1) + 1;

vector<vector<Ciphertext<DCRTPoly>>> EncryptedConstDivCoefficients;
vector<vector<Ciphertext<DCRTPoly>>> EncryptedConstEqCoefficients;

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
	vector<vector<ll>> ConstDivCoefficients(p, vector<ll>(p));
	
	for(int d = 1; d < p; ++d) {
		vector<ll> data_x(p);
		vector<ll> data_y(p);

		for(int x = 0; x < p; ++x) {
			data_x[x] = x;
			data_y[x] = x / d;
		}

		ConstDivCoefficients[d] = lagrange_interpolation(data_x, data_y);
	}

	for(int d = 0; d < p; ++d) {
		vector<Ciphertext<DCRTPoly>> EncryptedConstDivCoefficients_d;
		for(int x = 0; x < p; ++x) {
				Plaintext p_temp = cryptoContext->MakePackedPlaintext(vector<int64_t>(1, ConstDivCoefficients[d][x]));
    			Ciphertext<DCRTPoly> c_temp = cryptoContext->Encrypt(keyPair.publicKey, p_temp);
				EncryptedConstDivCoefficients_d.push_back(c_temp);
			}
		EncryptedConstDivCoefficients.push_back(EncryptedConstDivCoefficients_d);
	}
}

// g(x) := 1 :if(x=y), 0 :if(x!=y) とし、f(x) = g(x)を満たす関数f(x)について多項式補完を行い、その係数を求める
// d in Z_p について計算する。
void PrecomupteConstEqCoefficient() {
	vector<vector<ll>> ConstEqCoefficients(p, vector<ll>(p));
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

	for(int d = 0; d < p; ++d) {
		vector<Ciphertext<DCRTPoly>> EncryptedConstEqCoefficients_d;
		for(int x = 0; x < p; ++x) {
				Plaintext p_temp = cryptoContext->MakePackedPlaintext(vector<int64_t>(1, ConstEqCoefficients[d][x]));
    			Ciphertext<DCRTPoly> c_temp = cryptoContext->Encrypt(keyPair.publicKey, p_temp);
				EncryptedConstEqCoefficients_d.push_back(c_temp);
			}
		EncryptedConstEqCoefficients.push_back(EncryptedConstEqCoefficients_d);
	}
}


// Pows関数ではべき乗列の作成を工夫することで、乗算が深くなることを防いでいる
vector<Ciphertext<DCRTPoly>> Pows(const Ciphertext<DCRTPoly>& c) {
	int l = log2(p);
	vector<Ciphertext<DCRTPoly>> c_pows;

	// c_powsをp - 1個分の配列にresizeする
	for (int i = 0; i < p - 1; ++i) c_pows.push_back(c);

	// c^1からc^l-1までを計算する
	for (int i = 0; i < l; ++i) {
		int base = pow(2, i);

		for (int j = 0; j < base; ++j) {
			c_pows[base + j] = cryptoContext->EvalMult(c_pows[base - 1], c_pows[j]);
		}
	}

	// c^lからc^p-1までを計算する
	if (pow(2, l) < p - 1) {
		long long base = pow(2, l);

		for (int i = 0; i < p - base - 1; ++i) {
			c_pows[base + i] = cryptoContext->EvalMult(c_pows[base - 1], c_pows[i]);
		}
	}

	return c_pows;
}

Ciphertext<DCRTPoly> ConstDiv(vector<Ciphertext<DCRTPoly>> C_pow, int divisor) {
	// C_div_result = C_0 を作成する
	Plaintext p_div_result = cryptoContext->MakePackedPlaintext(vector<int64_t>(1, 0));
    Ciphertext<DCRTPoly> C_div_result = cryptoContext->Encrypt(keyPair.publicKey, p_div_result);

	for(int i = 0; i < p; ++i) {
		// x^0 = 1のため、定数のみ加算する
		if(i == 0) {
			C_div_result = cryptoContext->EvalAdd(C_div_result, EncryptedConstDivCoefficients[divisor][i]);
			continue;
		}

		Ciphertext<DCRTPoly> temp = cryptoContext->EvalMult(C_pow[i - 1], EncryptedConstDivCoefficients[divisor][i]);
		C_div_result = cryptoContext->EvalAdd(C_div_result, temp);
	}

	return C_div_result;
}

Ciphertext<DCRTPoly> ConstEq(vector<Ciphertext<DCRTPoly>> C_pow, int divisor) {
	// C_div_result = C_0 を作成する
	Plaintext p_div_result = cryptoContext->MakePackedPlaintext(vector<int64_t>(1, 0));
    Ciphertext<DCRTPoly> C_div_result = cryptoContext->Encrypt(keyPair.publicKey, p_div_result);

	for(int i = 0; i < p; ++i) {
		// x^0 = 1のため、定数のみ加算する
		if(i == 0) {
			C_div_result = cryptoContext->EvalAdd(C_div_result, EncryptedConstEqCoefficients[divisor][i]);
			continue;
		}

		Ciphertext<DCRTPoly> temp = cryptoContext->EvalMult(C_pow[i - 1], EncryptedConstEqCoefficients[divisor][i]);
		C_div_result = cryptoContext->EvalAdd(C_div_result, temp);
	}

	return C_div_result;
}

Ciphertext<DCRTPoly> Div(Ciphertext<DCRTPoly> C_a, Ciphertext<DCRTPoly> C_d) {
	// C_sum = C_0 を作成する
	Plaintext p_sum = cryptoContext->MakePackedPlaintext(vector<int64_t>(1, 0));
    Ciphertext<DCRTPoly> C_sum = cryptoContext->Encrypt(keyPair.publicKey, p_sum);

	// 各暗号文のべき乗列を計算する
	vector<Ciphertext<DCRTPoly>> C_a_pow = Pows(C_a);
	vector<Ciphertext<DCRTPoly>> C_d_pow = Pows(C_d);
	
	for(int i = 1; i < p; ++i) {
		// C_a/i の準同型演算の結果を返す
		Ciphertext<DCRTPoly> C_div_a_i = ConstDiv(C_a_pow, i);
		// d == i の場合：C_0
		// d != i の場合：C_1
		Ciphertext<DCRTPoly> C_equal_d_i = ConstEq(C_d_pow, i);

		// C_sum = C_sum + FHE.Mult(C_a/i, C_d=i)
		// d == i の場合：C_sum = C_sum + C_a/i * C_0 = C_sum + C_0  = C_sum (C_0より0を足しているだけなので、C_sumのまま)
		// d != i の場合：C_sum = C_sum + C_a/i * C_1 = C_sum + C_a/i        (C_a/iの結果がC_sumに加算される)
		Ciphertext<DCRTPoly> temp = cryptoContext->EvalMult(C_div_a_i, C_equal_d_i); // error zone
		C_sum = cryptoContext->EvalAdd(C_sum, temp);
	}

	Ciphertext<DCRTPoly> C_div_a_d = C_sum;

    return C_div_a_d;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
		cout << "few argument error" << endl;
		return 1;
	}
    
    // 各種パラメータの設定
    CCParams<CryptoContextBGVRNS> parameters;
    parameters.SetMultiplicativeDepth(multiplicative_depth);
    parameters.SetPlaintextModulus(p);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetStandardDeviation(3.2);
    parameters.SetKeySwitchTechnique(BV);
    parameters.SetDigitSize(1);
    parameters.SetRingDim(8);
    parameters.SetScalingTechnique(FIXEDAUTO);

    // 使用する機能の設定
    cryptoContext = GenCryptoContext(parameters);
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);

    // 鍵
    keyPair = cryptoContext->KeyGen(); // 秘密鍵・公開鍵を作成する。
    cryptoContext->EvalMultKeyGen(keyPair.secretKey); // 準同型乗算を行う際に使用する鍵を作成する。
    
    // 平文
    int64_t a = *argv[1] - '0', d = *argv[2] - '0';
	Plaintext plaintext_a = cryptoContext->MakePackedPlaintext(vector<int64_t>(1, a));
	Plaintext plaintext_d = cryptoContext->MakePackedPlaintext(vector<int64_t>(1, d));

    // 暗号化
    Ciphertext<DCRTPoly> cihpertext_a = cryptoContext->Encrypt(keyPair.publicKey, plaintext_a);
    Ciphertext<DCRTPoly> cihpertext_d = cryptoContext->Encrypt(keyPair.publicKey, plaintext_d);

    PrecomupteConstDivCoefficient();
	PrecomupteConstEqCoefficient();

	// 時間計測
	chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

	Ciphertext<DCRTPoly> C_div_result = Div(cihpertext_a, cihpertext_d);
	
	end = chrono::system_clock::now();
	double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);

	Plaintext div_result;
    cryptoContext->Decrypt(keyPair.secretKey, C_div_result, &div_result);

    cout << div_result << endl;

	printf("time %lf[ms]\n", time);
}