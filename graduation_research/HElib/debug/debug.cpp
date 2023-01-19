#include <iostream>
#include <cmath>
#include <helib/helib.h>
using namespace std;
using namespace helib;


// パラメータ関係
const int B_range = 4; // [-B^B_range, B^B_range]
const int ckks_m = 65536;// 131072
const int ckks_bits = 613;// 1445
const int ckks_precision = 40;// 40
const int ckks_c = 3; // 8

// 多項式で使用するデータ点
const int N = 9;

// 多項式補間なので、最大次数は データ点　- 1
const int P = N - 1;

// 逆関数 y = 1/x の近似範囲[1, A]
const int A = 10;

// B: 指数関数の入力値の範囲 [-B, B]
// B: 値は2のべき乗にすること(繰り返し二乗法で効率化するため)
// この値を大きくしすぎると、繰り返し二乗法の途中で乗算の深さの制限を飛び出してしまい値がめちゃくちゃになります。
const int B = pow(2, B_range);

// 指数関数 y = exp(x) の範囲[-1, 1]をラグランジュ補間した多項式の係数
// 左からx^0, x^1, ..., x^p と続く
// 係数値が小さいため、小数点3桁分を使用し残りは切り捨てる
//const double exp_lagrange_coefficients[] = {1, 0.999,  0.500, 0.167, 0.0505, 0.0113, -0.0312, -0.00325, 0.0227};

// 参考：全桁
const double exp_lagrange_coefficients[] = {1, 0.99949524,  0.50095556, 0.16737778, 0.05057778, 0.01137778, -0.03128889, -0.00325079, 0.02275556};


// 逆関数 y = 1/x の範囲[1, a]をラグランジュ補間した係数
// 左からx^0, x^1, ..., x^p と続く
// 小数点3桁までを使用し残りは切り捨てる
//const double inverse_lagrange_coefficients[] = {2.659, -2.811, 1.589, -0.538, 0.114, -1.549e-2, 1.287e-3, -5.995e-5, 1.198e-6};

// 参考：全桁
const double inverse_lagrange_coefficients[] = {2.65967021,   -2.81125582,     1.58988512,    -0.538805990, 0.114775235, -1.54970593e-2, 1.28706134e-3, -5.99599246e-5,  1.19861852e-6};

// input : c
// output: c_pows = { c^1, c^2, ..., c^{p-1} } 
// Pows関数ではべき乗列の算出を工夫することで、乗算が深くなることを防いでいる
vector<helib::Ctxt> Pows(helib::Ctxt c) {
	int l = log2(P);
	vector<helib::Ctxt> c_pows;
	
	// c_powsをp個分の配列にresizeする
	for (int i = 0; i < P; ++i) c_pows.push_back(c);

	// c^1からc^l-1までを計算する
	for (int i = 0; i < l; ++i) {
		int base = pow(2, i);
		for (int j = 0; j < base; ++j) {
			c_pows[base + j] = c_pows[base - 1];
			c_pows[base + j] *= c_pows[j];
		}
	}

	// c^lからc^p-1までを計算する
	if (pow(2, l) < P - 1) {
		long long base = pow(2, l);

		for (int i = 0; i < P - base - 1; ++i) {
			c_pows[base + i] = c_pows[base - 1];
			c_pows[base + i] *= c_pows[i];
		}
	}

	return c_pows;
}


// y = exp(x/B) ** B
// x: 指数関数の引数
helib::Ctxt HE_exp(helib::Ctxt x, const helib::Context& context, const helib::PubKey& pk) {
	// C_exp_x = C_0 を作成する
	vector<double> init = {0};
	PtxtArray P_exp_x(context, init);
	Ctxt C_exp_x(pk);
	P_exp_x.encrypt(C_exp_x);

	// スケーリング
	x *= 1.0/B;

	vector<helib::Ctxt> C_pows = Pows(x);

	// ラグランジュ補間で求めた多項式を算出する
	for(int i = 0; i < N; ++i) {
		// x^0では係数を足すだけ
		if(i == 0) {
			C_exp_x += exp_lagrange_coefficients[i];

			// C_inverse_exp_xの確認
			/*
			PtxtArray P_inverse_exp_x(context);
			P_inverse_exp_x.decrypt(C_inverse_x, sk);
			vector<double> sigmoid_vaue;
			P_inverse_exp_x.store(sigmoid_vaue);
			cout << "x^" << i + 1 << ": " << sigmoid_vaue[0] << endl;
			*/
			continue;
		}

		// C_powsはC_pows[0] = C^1, C_pows[1] = C^2と、添え字とべき数がズレているため、i - 1 している
		Ctxt temp(C_pows[i - 1]);
		temp *= exp_lagrange_coefficients[i];
		C_exp_x += temp;

		// C_inverse_exp_xの確認
		/*
		PtxtArray P_inverse_exp_x(context);
		P_inverse_exp_x.decrypt(C_inverse_x, sk);
		vector<double> sigmoid_vaue;
		P_inverse_exp_x.store(sigmoid_vaue);
		cout << "x^" << i << ": " << sigmoid_vaue[0] << endl;
		*/
	}
	
	// 繰り返し二乗法で、y = exp(x/B)^B のB乗の部分を算出する
	for(int i = 0; i < floor(log2(B)); ++i) {
		C_exp_x *= C_exp_x;
	}
	
	return C_exp_x;
}

// y = 1/x [1 <= x] 1以上の値の逆数を求められる。0 ～1 未満はダメ
helib::Ctxt HE_inverse(helib::Ctxt x, const helib::Context& context, const helib::PubKey& pk, const helib::SecKey& sk) {
	// C_inverse_x = C_0 を作成する
	vector<double> init = {0};
	PtxtArray P_inverse_x(context, init);
	Ctxt C_inverse_x(pk);
	P_inverse_x.encrypt(C_inverse_x);

	// 入力値xを復号する
	// スケーリング値(xを近似範囲である[1, a]の間の値に変換する)
    // ！ここで「切り捨て」と「対数関数」という準同型暗号では演算できない関数が出てくる！
    // 現状はこのスケーリング値の算出時のみ、xを復号する必要がある
	PtxtArray x_client(context);
	x_client.decrypt(x, sk);
	vector<double> x_decrypt;
	x_client.store(x_decrypt);

	// A = 10の場合、常用対数で計算できる
	int R = -1;
	if (A == 10) {
		R = pow(A, floor(log10(x_decrypt[0])));
	}
	else {
		// a != 10の場合、底の変換公式を使う
		// log_a(x) = log_e(x) / log_e(a) 
		R = pow(A, floor(log(x_decrypt[0]) / log(A)));
	}
	
	x *= (1.0/R);
	vector<helib::Ctxt> C_pows = Pows(x);

	// powsの確認
	/*
	for(int i = 0; i < P; ++i) {
		PtxtArray P_inverse_exp_x(context);
		P_inverse_exp_x.decrypt(C_pows[i], sk);
		vector<double> sigmoid_vaue;
		P_inverse_exp_x.store(sigmoid_vaue);
		cout << "x^" << i + 1 << ": " << sigmoid_vaue[0] << endl;
	}
	*/

	// ラグランジュ補間で求めた多項式を算出する
	for(int i = 0; i < N; ++i) {
		// x^0では係数を足すだけ
		if(i == 0) {
			C_inverse_x += inverse_lagrange_coefficients[i];
			continue;
		}

		// C_powsはC_pows[0] = C^1, C_pows[1] = C^2と、添え字とべき数がズレているため、i - 1 している
		Ctxt temp(C_pows[i - 1]);
		temp *= inverse_lagrange_coefficients[i];
		C_inverse_x += temp;
	}

	// スケーリングした分をもとに戻す
	C_inverse_x *= (1.0 / R);

	return C_inverse_x;
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		cout << "few argument error" << endl;
		return 1;
	}

	Context context = ContextBuilder<CKKS>()
	// m: cyclotomic index
	// mを増加させるほど...
	//   メリット　：セキュリティとslotが増加する
	//   デメリット：パフォーマンスの低下し、暗号文のサイズが増加する
	.m(ckks_m)
	
	// bits: "ciphertext modulus"のビット数  
	// bitsを増加させるほど...
	//   メリット　：乗算の深さが増加する
	//   デメリット：セキュリティが下がる
	// 深さとbitsの関係については「HElib/examples/tutorial/02_depth.cpp」を参照すること
	.bits(ckks_bits)

	// precision: 暗号化・復号時の精度
	// precisionを増加させるほど...
	//   メリット　：精度が向上する
	//   デメリット：乗算の深さが減少する
	// この値は40以下で使用することが推奨されている
	.precision(ckks_precision)

	// c: "key-switching matrices"の列数
	// cを増加させるほど...
	//   メリット　：セキュリティが少し上がる
	//   デメリット：パフォーマンスが低下し、"public key"のメモリサイズが増加する
	// この値は8以下で使用することが推奨されている
	.c(ckks_c)
	.build();

	cout << "securityLevel=" << context.securityLevel() << "\n";
    cout << "nslots=" << context.getNSlots() << "\n";

	SecKey secretKey(context);
	secretKey.GenSecKey();

	const PubKey& publicKey = secretKey;
	
	string input = argv[1]; // 引数の受け取り
	double input_value = stod(input);
	vector<double> x = {-input_value}; // y = exp(-x) なので引数にマイナスを掛ける
	PtxtArray P_x(context, x);
	Ctxt C_x(publicKey);
	P_x.encrypt(C_x);

	// 時間計測
	chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

	// シグモイド関数 y = 1 / (1 + exp(-x)) の算出
	// ！注意！：x は [-B, B]の範囲内のみ(ここから外れるとexpを計算できない)
	
	cout << "C_x.capacity=" << C_x.capacity() << " ";
    cout << "C_x.errorBound=" << C_x.errorBound() << "\n";
	
	// exp(-x)
	helib::Ctxt C_exp_x = HE_exp(C_x, context, publicKey);
	/*
	cout << "C_exp_x.capacity=" << C_exp_x.capacity() << " ";
    cout << "C_exp_x.errorBound=" << C_exp_x.errorBound() << "\n";
	*/
	// 1 + exp(-x)
	C_exp_x += 1.0;

	// print 1 + exp(-x)
	/*
	PtxtArray pp(context);
	pp.decrypt(C_exp_x, secretKey);
	vector<double> vv;
	pp.store(vv);
	cout << "C_exp_x + 1 value=" << vv[0] << endl;
	*/

	
	cout << "C_exp_x + 1.capacity=" << C_exp_x.capacity() << " ";
    cout << "C_exp_x + 1.errorBound=" << C_exp_x.errorBound() << "\n";
	

	// C_exp_xを一度復号してノイズを取り除き、その後再び暗号化する
	
	PtxtArray P_re_exp_x(context);
	P_re_exp_x.decrypt(C_exp_x, secretKey);
	Ctxt C_re_exp_x(publicKey);
	P_re_exp_x.encrypt(C_re_exp_x);
	
	// 1 / (1 + exp(-x))
	helib::Ctxt C_inverse_exp_x = HE_inverse(C_re_exp_x, context, publicKey, secretKey);
	//helib::Ctxt C_inverse_exp_x = HE_inverse(C_exp_x, context, publicKey, secretKey);
	
	cout << "C_inverse_exp_x.capacity=" << C_inverse_exp_x.capacity() << " ";
    cout << "C_inverse_exp_x.errorBound=" << C_inverse_exp_x.errorBound() << "\n";
	
	end = chrono::system_clock::now();
	double time_measurement = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);

	PtxtArray P_inverse_exp_x(context);
	P_inverse_exp_x.decrypt(C_inverse_exp_x, secretKey);
	vector<double> sigmoid_vaue;
	P_inverse_exp_x.store(sigmoid_vaue);
	cout << input_value << ": " << sigmoid_vaue[0] << endl;

	printf("time %lf[ms]\n", time_measurement);

	/*
	// sigmoid の確認
	for(double i = -1; i <= 1; i += 0.1) {
		vector<double> x = {-i};
		PtxtArray P_x(context, x);
		Ctxt C_x(publicKey);
		P_x.encrypt(C_x);

		// exp(-x)
		helib::Ctxt C_exp_x = HE_exp(C_x, context, publicKey);

		// 1 + exp(-x)
		C_exp_x += 1.0;

		// 1 / (1 + exp(-x))
		helib::Ctxt C_inverse_exp_x = HE_inverse(C_exp_x, context, publicKey, secretKey);

		
	}
	*/

	/*
	// expの確認
	for(double i = -2; i <= 2; i += 0.5) {
		vector<double> x = {-i};
		PtxtArray P_x(context, x);
		Ctxt C_x(publicKey);
		P_x.encrypt(C_x);

		// exp(-x)
		helib::Ctxt C_exp_x = HE_exp(C_x, context, publicKey);
		PtxtArray pp(context);
		pp.decrypt(C_exp_x, secretKey);
		vector<double> vv;
		pp.store(vv);
		cout << i << ": " << vv[0] << endl;
	}
	*/

	/* 
	// HE_inverseの確認
	double test = 1.231;
	vector<double> x = {test};
	PtxtArray P_x(context, x);
	Ctxt C_x(publicKey);
	P_x.encrypt(C_x);
	helib::Ctxt result = HE_inverse(C_x, context, publicKey, secretKey);

	PtxtArray P_test(context);
	P_test.decrypt(result, secretKey);
	vector<double> t;
	P_test.store(t);
	cout << test << ": " << t[0] << endl;
	*/
}
