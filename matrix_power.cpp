#include <bits/stdc++.h>
#define int long long

using namespace std;

long long mod = 1e9 + 7;
typedef vector<int> row;
typedef vector<row> Matrix;
vector<Matrix> powers(64);

/// @brief Multiply two matrices a, and b
/// @param a
/// @param b
/// @return result of multiplication
Matrix multiply(Matrix &a, Matrix &b) {
    int n = a.size(), m = b[0].size();
    Matrix ret(n, row(m, 0));

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < a[0].size(); k++) {
            if (a[i][k] == 0)
                continue;

            for (int j = 0; j < m; j++)
                ret[i][j] = (ret[i][j] + a[i][k] * b[k][j] % mod) % mod;
        }
    }
    return ret;
}

/// @brief initialize all powers of a matrix
/// @param ret
void init(Matrix &ret) {
    powers[0] = ret;
    for (int i = 1; i < 64; i++)
        powers[i] = multiply(powers[i - 1], powers[i - 1]);
}

/// @brief raise a matrix to a power p
/// @param p
/// @return result of power
Matrix get_power(int p) {
    int n = powers[0].size();
    Matrix ret(n, row(n, 0));
    for (int i = 0; i < n; i++)
        ret[i][i] = 1;

    for (int i = 0; i < 64; i++) {
        if ((p >> i) & 1)
            ret = multiply(ret, powers[i]);
    }

    return ret;
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
