#include <bits/stdc++.h>

using namespace std;
using ll = long long;

const int K = 60, N = 505;

struct xor_basis {
    int sz;
    ll basis[K];
    bitset<N> index[K];

    xor_basis(): sz(0) {
        for (auto& i : index)
            i.reset();
        memset(basis, 0, sizeof basis);
    }

    void insert(ll mask, int idx) {
        bitset<N> cur; cur[idx] = 1;
        for (int i = K - 1; i >= 0; i--) {
            if (!(mask >> i & 1))
                continue;
            if (!basis[i]) {
                basis[i] = mask, sz++;
                index[i] = cur;
                return;
            }
            mask ^= basis[i];
            cur ^= index[i];
        }
    }
    bool check(ll mask) {
        for (int i = K - 1; i >= 0; i--) {
            if (!(mask >> i & 1))
                continue;
            if (!basis[i]) {
                return false;
            }
            mask ^= basis[i];
        }
        return true;
    }
    ll calc_xors() {
        return (1ll << sz);
    }
    ll max() {
        ll ans = 0;
        for (int i = K - 1; i >= 0; i--) {
            if (!basis[i])continue;
            if (ans >> i & 1)continue;
            ans ^= basis[i];
        }
        return ans;
    }
    ll kth(ll k) {
        ll tot = (1ll << sz), ans = 0;
        for (int i = K - 1; i >= 0; i--) {
            if (!basis[i])continue;

            tot /= 2;
            if (tot < k && !(ans >> i & 1))
                ans ^= basis[i];
            else if (tot >= k && (ans >> i & 1))
                ans ^= basis[i];

            if (k > tot)
                k -= tot;
        }
        return ans;
    }
    bool construct(ll x, bitset<N>& ret) {
        for (int i = K - 1; i >= 0; i--) {
            if (!(x >> i & 1))
                continue;
            if (!basis[i])
                return false;

            x ^= basis[i];
            ret ^= index[i];
        }
        return x == 0;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n = 0, q;
    cin >> q;

    xor_basis b;

    while (q--) {
        int t; ll x;
        cin >> t;

        if (t == 1) { // insert x
            cin >> x;
            b.insert(x, n++);
        }
        else if (t == 2) { // find x-th distinct xor
            cin >> x;
            cout << b.kth(x) << "\n";
        }
        else if (t == 3) { // calc number of ways to get x as xor
            cin >> x;
            if (!b.check(x)) {
                cout << 0 << "\n";
                continue;
            }
            ll ans = (1ll << (n - b.sz));
            cout << ans << "\n";
        }
        else if (t == 4) { // calc max xor value
            cout << b.max() << "\n";
        }
        else if (t == 5) { // construct a sub-sequence with xor = x
            cin >> x;
            bitset<N> ret; ret.reset();
            bool ok = b.check(x);
            if (!ok) {
                cout << -1 << "\n";
                continue;
            }
            for (int i = 0; i < n; i++)
                cout << ret[i];
            cout << "\n";
        }
    }

    return 0;
}
