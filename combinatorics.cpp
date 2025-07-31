#include <bits/stdc++.h>

using namespace std;

#define int long long
const int N = 1e6 + 20, mod = 1e9 + 7;
int fact[N], ifact[N];

int add(int a, int b) {
    a += b;
    a %= mod;
    if (a < 0)
        a += mod;
    return a;
}

int mul(int a, int b) {
    a *= b;
    a %= mod;
    return a;
}

int fpow(int a, int b) {
    if (b == 0)
        return 1;
    int ret = fpow(a, b / 2);
    ret = mul(ret, ret);
    if (b & 1)
        ret = mul(ret, a);
    return ret;
}

int inv(int a) {
    return fpow(a, mod - 2);
}

int nCk(int n, int k) {
    return mul(fact[n], mul(ifact[k], ifact[n - k]));
}

void pre() {
    memset(fact, 1, sizeof fact);
    for (int i = 1; i < N; i++)
        fact[i] = mul(i, fact[i - 1]);
    for (int i = 1; i < N; i++)
        ifact[i] = inv(fact[i]);
    
    ifact[N - 1] = inv(fact[N - 1]);
    for (int i = N - 2; i >= 0; i--)
        ifact[i] = mul(i + 1, ifact[i + 1]);
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
