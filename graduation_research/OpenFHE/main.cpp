#include <openfhe.h>
using namespace std;
using namespace lbcrypto;

int main() {
    // multicative depth
    uint32_t multDepth = 4;

    // scale size
    // 実数値に係数を掛けて整数に近似する際の、"係数"を決めるのがこのパラメータ
    // この値は精度とノイズのトレードオフとなる。
    // 係数を大きくするほど精度は上がるが、演算によるノイズ量も増えてしまう
    uint32_t scaleModSize = 50;

    // slot size
    // slotの数は GetRingDimensin() の値に依存している
    // この値に関係なく指定することもできる
    uint32_t batchSize = 1;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleModSize);
    parameters.SetBatchSize(batchSize);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    std::vector<double> x1 = {3};

    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x1);

    auto c1 = cc->Encrypt(keys.publicKey, ptxt1);

    vector<Ciphertext<DCRTPoly>> c_pows;

    int n = 6;
    // c_powsをp - 1個分の配列にresizeする
    for (int i = 0; i < n; ++i) c_pows.push_back(c1);
    for (int i = 1; i < n; ++i) c_pows[i] = cc->EvalMult(c_pows[i], c_pows[i - 1]);
    for (int i = 0; i < n; ++i) {
        Plaintext result;
        std::cout.precision(8);
        cc->Decrypt(keys.secretKey, c_pows[i], &result);
        cout << cc->MakePlaintext() << endl;
    }

    /*
    std::vector<double> x1 = {-10.5};
    std::vector<double> x2 = {2.1};

    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x1);
    Plaintext ptxt2 = cc->MakeCKKSPackedPlaintext(x2);

    auto c1 = cc->Encrypt(keys.publicKey, ptxt1);
    auto c2 = cc->Encrypt(keys.publicKey, ptxt2);

    auto cMul = cc->EvalMult(c1, c2);
    auto abc = cc->EvalMult(c1, 2.5);
    Plaintext result;
    std::cout.precision(8);

    cc->Decrypt(keys.secretKey, c1, &result);
    cout << result << endl;

    cc->Decrypt(keys.secretKey, c2, &result);
    cout << result << endl;

    cc->Decrypt(keys.secretKey, cMul, &result);
    cout << result << endl;

    cc->Decrypt(keys.secretKey, abc, &result);
    cout << result << endl;
    */
}
