#include <bits/stdc++.h>

using namespace std;
using ll = long long;

// Hashing
mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
int rand(int l, int r) { return uniform_int_distribution<int>(l, r)(rng); }

struct Hash {
    int mod, base, n;
    vector<int> prf, suf, pows, invs;
    vector<int> primesMod = {
        1000000007, 1000000009, 1000000021, 1000000033, 1000000087, 1000000093, 1000000097, 1000000103, 1000000123,
        1000000181, 1000000207, 1000000223, 1000000241, 1000000271, 1000000289, 1000000297, 1000000321, 1000000349,
        1000000363, 1000000403, 1000000409, 1000000411, 1000000427, 1000000433, 1000000439, 1000000447, 1000000453,
        1000000459, 1000000483, 1000000513, 1000000531, 1000000579, 1000000607, 1000000613, 1000000637, 1000000663,
        1000000711, 1000000753, 1000000787, 1000000801, 1000000829, 1000000861, 1000000871, 1000000891, 1000000901,
        1000000919, 1000000931, 1000000933, 1000000993, 1000001011
    };

    int fpow(int x, int p) {
        int ans = 1;
        while (p) {
            if (p & 1)
                ans = (1LL * ans * x) % mod;
            p /= 2;
            x = (1LL * x * x) % mod;
        }
        return ans;
    }

    Hash(const string &s) {
        int r = rand(0, 49);
        mod = primesMod[r];
        base = rand(2, 199);

        n = s.size();
        pows.assign(n + 5, 1);
        invs.assign(n + 5, 1);

        for (int i = 1; i < n + 5; i++) {
            pows[i] = (1LL * pows[i - 1] * base) % mod;
            invs[i] = fpow(pows[i], mod - 2);
        }

        int value = 0, p = 0;
        prf.push_back(0);
        for (int i = 0; i < n; i++) {
            char c = s[i];
            value = (value + 1LL * (int)c * pows[p++]) % mod;
            prf.push_back(value);
        }
        prf.push_back(0);

        value = 0, p = 0;
        suf.push_back(0);
        for (int i = n - 1; i >= 0; i--) {
            char c = s[i];
            value = (value + 1LL * (int)c * pows[p++]) % mod;
            suf.push_back(value);
        }
    }

    int substring(int l, int r, int who = 0) /// 1-based
    {
        int cand;
        if (who == 0) /// prfix
        {
            cand = (prf[r] - prf[l - 1]) % mod;
            cand = (cand + mod) % mod;
            cand = (1LL * cand * invs[l]) % mod;
        } else /// suffix
        {
            int ll = l, rr = r;
            l = n - rr + 1, r = n - ll + 1;
            cand = (suf[r] - suf[l - 1]) % mod;
            cand = (cand + mod) % mod;
            cand = (1LL * cand * invs[l]) % mod;
        }
        return cand;
    }

    bool palindrome(int l, int r) {
        return substring(l, r, 0) == substring(l, r, 1);
    }
};

// ==========================================================================================

// 2D Hashing
// mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
// int rand(int l, int r) { return uniform_int_distribution<int>(l, r)(rng); }

vector<int> primesMod = {
    1000000007, 1000000009, 1000000021, 1000000033, 1000000087,
    1000000093, 1000000097, 1000000103, 1000000123, 1000000181,
    1000000207, 1000000223, 1000000241, 1000000271, 1000000289,
    1000000297, 1000000321, 1000000349, 1000000363, 1000000403,
    1000000409, 1000000411, 1000000427, 1000000433, 1000000439,
    1000000447, 1000000453, 1000000459, 1000000483, 1000000513,
    1000000531, 1000000579, 1000000607, 1000000613, 1000000637,
    1000000663, 1000000711, 1000000753, 1000000787, 1000000801,
    1000000829, 1000000861, 1000000871, 1000000891, 1000000901,
    1000000919, 1000000931, 1000000933, 1000000993, 1000001011
};
int fpow(int x, int p) {
    int ans = 1;
    while (p) {
        if (p & 1)
            ans = (1LL * ans * x) % mod;
        p /= 2;
        x = (1LL * x * x) % mod;
    }
    return ans;
}

struct Hash2d {
    int n, m;
    int r, mod, X, Y;
    vector<vector<int> > A;
    vector<int> px, ipx, py, ipy;

    Hash2d(const vector<string> &a) {
        r = rand(0, 49);
        mod = primesMod[r];
        X = rand(2, 100);
        Y = rand(101, 199);

        n = a.size(), m = a[0].size();
        px.resize(n + 1, 1);
        ipx.resize(n + 1, 1);
        py.resize(m + 1, 1);
        ipy.resize(m + 1, 1);

        for (int i = 1; i <= n; i++) {
            px[i] = (px[i - 1] * X) % mod;
            ipx[i] = fpow(px[i], mod - 2);
        }
        for (int i = 1; i <= m; i++) {
            py[i] = (py[i - 1] * Y) % mod;
            ipy[i] = fpow(py[i], mod - 2);
        }

        A = vector<vector<int> >(n + 1, vector<int>(m + 1, 0));
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++) {
                A[i][j] = (int)a[i - 1][j - 1];
                A[i][j] = (A[i][j] * px[i]) % mod;
                A[i][j] = (A[i][j] * py[j]) % mod;
            }
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++)
                A[i][j] = (A[i][j] + A[i][j - 1]) % mod;
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++)
                A[i][j] = (A[i][j] + A[i - 1][j]) % mod;
    }
    int subgrid(int i1, int j1, int i2, int j2) {
        assert(1 <= i1 && i1 <= i2 && i2 <= n);
        assert(1 <= j1 && j1 <= j2 && j2 <= m);

        int sum = 0;
        sum += A[i2][j2];
        sum -= A[i1 - 1][j2];
        sum -= A[i2][j1 - 1];
        sum += A[i1 - 1][j1 - 1];
        sum %= mod;
        sum = (sum + mod) % mod;

        sum = (sum * ipx[i1 - 1]) % mod;
        sum = (sum * ipy[j1 - 1]) % mod;

        return sum;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    return 0;
}
