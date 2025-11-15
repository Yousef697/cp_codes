#include <bits/stdc++.h>

using namespace std;
using ll = long long;

struct basis {
    int n = 0;
    vector<ll> base;

    basis() {}
    basis(const vector<ll>& v) {
        for (auto i : v) {
            add(i);
        }
    }

    void add(ll x) {
        ll a = x;
        for (auto j : base) {
            a = min(a, a ^ j);
        }
        if (a) {
            n++;
            base.push_back(a);
        }
    }
    bool can(ll x) {
        for (auto i : base)
            x = min(x, x ^ i);
        return x == 0;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n;
    cin >> n;

    vector<ll> a(n);
    for (auto& i : a)
        cin >> i;

    basis b(a);
    cout << (1ll << b.n) << "\n";
    
    return 0;
}
