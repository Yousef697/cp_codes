#include <bits/stdc++.h>
#define int long long

using namespace std;

long long mod = 1e9 + 7;
typedef vector<int> row;
typedef vector<row> Matrix;

Matrix multiply(Matrix &a, Matrix &b) {
    int n = a.size(), m = b[0].size();
    Matrix ret(n, row(m, 0));

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < a[0].size(); k++) {
            if (a[i][k] == 0)
                continue;

            for (int j = 0; j < m; j++) {
                ret[i][j] = (ret[i][j] + a[i][k] * b[k][j] % mod);
                if (ret[i][j] >= mod)
                    ret[i][j] -= mod;
            }
        }
    }
    return ret;
}

Matrix get_power(Matrix& a, int p) {
    int n = a.size();
    Matrix ret(n, row(n, 0));
    for (int i = 0; i < n; i++)
        ret[i][i] = 1;

    while (p) {
        if (p & 1)
            ret = multiply(ret, a);
        a = multiply(a, a);
        p >>= 1;
    }

    return ret;
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
