#include <bits/stdc++.h>
#define int long long

using namespace std;

vector<int> add(vector<int> a, vector<int> b) {
    int n = a.size();
    for (int i = 0; i < n; i++)
        a[i] += b[i];
    return a;
}

vector<int> subtract(vector<int> a, vector<int> b) {
    int n = a.size();
    for (int i = 0; i < n; i++)
        a[i] -= b[i];
    return a;
}

vector<int> karatsuba(const vector<int>& p, const vector<int>& q) {
    int n = p.size();
    if (n <= 2) {
        vector<int> ret(2 * n, 0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                ret[i + j] += p[i] * q[j];
        return ret;
    }

    vector<int> pl(n / 2), pr(n / 2);
    vector<int> ql(n / 2), qr(n / 2);
    for (int i = 0; i < n / 2; i++) {
        pl[i] = p[i], ql[i] = q[i];
        pr[i] = p[i + n / 2], qr[i] = q[i + n / 2];
    }

    vector<int> z1 = karatsuba(pl, ql);
    vector<int> z2 = karatsuba(pr, qr);
    vector<int> z3 = karatsuba(add(pl, pr), add(ql, qr));
    z3 = subtract(subtract(z3, z1), z2);

    vector<int> ans(2 * n);
    for (int i = 0; i < n; i++) {
        ans[i] += z1[i];
        ans[i + n] += z2[i];
        ans[i + n / 2] += z3[i];
    }
    return ans;
}

vector<int> multiply(vector<int> p, vector<int> q) {

    int n = 1;
    while (n < max(p.size(), q.size()))
        n *= 2;

    while (p.size() < n)
        p.push_back(0);
    while (q.size() < n)
        q.push_back(0);

    vector<int> ret = karatsuba(p, q);
    return ret;
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n, m;

    cin >> n;
    vector<int> a(n);
    for (auto& i : a)
        cin >> i;
    
    cin >> m;
    vector<int> b(m);
    for (auto& j : b)
        cin >> j;

    vector<int> ret = multiply(a, b);
    for (int i = 0; i < n + m - 1; i++)
        cout << ret[i] << " ";

    return 0;
}
