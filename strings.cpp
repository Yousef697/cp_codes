#include <bits/stdc++.h>
#define int long long

using namespace std;

int mod = 1e9 + 7;

// KMP, String Matching
/*
    prefix: any sub string that starts from index 0
    suffix: any sub string that ends from index size-1
    proper prefix: any prefix that its size is not equal to the main string size

    compute_prefixes function (failure function) compute the proper prefixes of a string such that
    ret[i] = max(proper prefix: prroper prefix = suffix )
*/
vector<int> compute_prefixes(const string& pat) {

    int m = pat.size(), i, k;
    vector<int> prefixes(m);

    for (i = 1, k = 0; i < m; i++) {
        while (k > 0 && pat[k] != pat[i])
            k = prefixes[k - 1]; // return to the promising position

        if (pat[i] == pat[k])
            prefixes[i] = ++k;
        else
            prefixes[i] = k;
    }

    return prefixes;
}
vector<vector<int>> compute_automation(string s) {
    s += '#';
    int n = s.size();
    auto pi = compute_prefixes(s);
    vector<vector<int>> automation(n, vector<int>(26, 0));

    for (int i = 0; i < n; i++) {
        for (int ch = 0; ch < 26; ch++) {
            if (i > 0 && 'a' + ch != s[i])
                automation[i][ch] = automation[pi[i - 1]][ch];
            else
                automation[i][ch] = i + ('a' + ch == s[i]);
        }
    }
    return automation;
}
int KMP(const string& s, const string& pat) {

    int n = s.size(), m = pat.size(), ans = 0;
    vector<int> prefixes = compute_prefixes(pat);

    int i; // for s
    int k; // for pat

    for (i = 0, k = 0; i < n; i++) {
        while (k > 0 && pat[k] != s[i])
            k = prefixes[k - 1]; // back to the promising position

        if (s[i] == pat[k])
            k++;
        if (k == m)
            ans++, k = prefixes[k - 1]; // make k fail
    }

    return ans;
}

// ==========================================================================================

// Z-Algorithm
vector<int> z_algorithm(const string &s) {
    int n = s.size();
    vector<int> z(n);
    int l = 0, r = 0;

    for (int i = 1; i < n; i++) {
        if (i < r)
            z[i] = min(r - i, z[i - l]);

        while (i + z[i] < n && s[i + z[i]] == s[z[i]])
            z[i]++;

        if (i + z[i] > r)
            l = i, r = i + z[i];
    }

    return z;
}
vector<int> find_pattern(const string &s, const string &t) {
    string A = t + '#' + s;
    int n = s.size(), m = t.size();
    vector<int> Z = z_algorithm(A);

    vector<int> ans;
    for (int i = m + 1; i < n + m + 1; i++) {
        if (Z[i] == m)
            ans.push_back(i);
    }

    return ans;
}
int distinct_substrings(const string &s) {
    int m = 0, ans = 0;
    string t = "";

    for (const char &c: s) {
        m++;
        t += c;
        int mx = 0;
        reverse(t.begin(), t.end());

        vector<int> Z = z_algorithm(t);
        for (const int i: Z)
            mx = max(mx, i);
        ans += m - mx;

        reverse(t.begin(), t.end());
    }

    return ans;
}
string period(const string &s) {
    int n = s.size(), ans = 0;
    vector<int> Z = z_algorithm(s);

    for (int i = 1; i <= n; i++) {
        if (n % i)
            continue;

        if (i + Z[i] == n) {
            ans = i;
            break;
        }
    }

    return s.substr(0, ans);
}

// ==========================================================================================

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

// ==========================================================================================

// Recursive Trie
const int K = 2;
struct Trie {
    vector<int> leaf;
    vector<int> cnts;
    vector<vector<int> > tree;

    Trie() {
        leaf.push_back(0);
        cnts.push_back(0);
        tree.push_back(vector<int>(K, 0));
    }

    void insert(const string &s, int v = 0, int i = 0) {
        if (i == s.size()) {
            cnts[v]++;
            leaf[v] = 1;
            return;
        }

        if (!tree[v][s[i] - '0']) {

            tree[v][s[i] - '0'] = tree.size();
            leaf.push_back(0);
            cnts.push_back(0);
            tree.push_back(vector<int>(2, 0));
        }

        int next = tree[v][s[i] - '0'];
        insert(s, next, i + 1);
        cnts[v]++;
    }
    void erase(const string &s, int v = 0, int i = 0) {
        if (i == s.size()) {
            cnts[v]--;
            if (cnts[v] == 0)
                leaf[v] = 0;
            return;
        }

        assert(tree[v][s[i] - '0']);
        int next = tree[v][s[i] - '0'];
        erase(s, next, i + 1);

        cnts[v]--;
        if (cnts[next] == 0)
            tree[v][s[i] - '0'] = 0;
    }
    int serach_word(const string &s, int v = 0, int i = 0) {
        if (i == s.size())
            return leaf[v];

        if (!tree[v][s[i] - '0'])
            return -1;
        return serach_word(s, tree[v][s[i] - '0'], i + 1);
    }
    int serach_prefix(const string &s, int v = 0, int i = 0) {
        if (i == s.size())
            return 1;

        if (!tree[v][s[i] - '0'])
            return -1;
        return serach_prefix(s, tree[v][s[i] - '0'], i + 1);
    }
};

// ==========================================================================================

// Manacher
vector<int> manacher(const string &s) {
    string t = "";
    for (const char &c: s)
        t += '#', t += c;
    t += '#';
    int n = t.size();
    t = '!' + t + '*';

    int l = 0, r = 1;
    vector<int> p(n + 2);
    for (int i = 1; i <= n; i++) {
        p[i] = min(r - i, p[l + (r - i)]);
        p[i] = max(0ll, p[i]);

        while (t[i - p[i]] == t[i + p[i]])
            p[i]++;

        if (i + p[i] > r)
            l = i - p[i], r = i + p[i];
    }

    for (int &i: p)
        i--;

    // [even, odd]
    return vector<int>(p.begin() + 1, p.end() - 1);
}

// ==========================================================================================

// Suffix Array
struct SuffixArray {

    int n;
    string s;
    vector<int> suf, classes, lcp;

    SuffixArray(string& str) {
        s = str + char(0), n = s.size();
        build();
        calc_lcp();
    }

    void count_sort(vector<int> &suf, vector<int> &classes) {
        int n = suf.size(), d = 2;
        vector<int> freq(n), pos(n);
        for (int i = 0; i < n; i++) {
            freq[classes[i]]++;
        }

        pos = freq;
        for (int i = 1; i < n; i++) {
            pos[i] += pos[i - 1];
        }

        vector<int> new_suf(n);
        for (int i = n - 1; i >= 0; i--) {
            int x = classes[suf[i]];
            pos[x]--;
            new_suf[pos[x]] = suf[i];
        }
        suf = new_suf;
    }
    void build() {
        suf.resize(n);
        classes.resize(n);

        // k = 0
        vector<pair<char, int> > a(n);
        for (int i = 0; i < n; i++) {
            a[i] = {s[i], i};
        }
        sort(a.begin(), a.end());

        for (int i = 0; i < n; i++) {
            suf[i] = a[i].second;
        }
        for (int i = 1; i < n; i++) {
            classes[suf[i]] = classes[suf[i - 1]] + (a[i].first != a[i - 1].first);
        }

        int k = 0;
        while ((1 << k) < n) {
            // k -> k + 1
            for (int i = 0; i < n; i++) {
                suf[i] = (suf[i] - (1 << k) + n) % n;
            }
            count_sort(suf, classes);

            vector<int> c(n);
            c[suf[0]] = 0;
            for (int i = 1; i < n; i++) {
                pair<int, int> cur = {classes[suf[i]], classes[(suf[i] + (1 << k)) % n]};
                pair<int, int> prv = {classes[suf[i - 1]], classes[(suf[i - 1] + (1 << k)) % n]};
                c[suf[i]] = c[suf[i - 1]] + (cur != prv);
            }
            classes = c;
            k++;
        }
    }
    int comp(int i, const string& t) {
        int j = 0;
        while (i < n && j < t.size() && s[i] == t[j]) {
            i++, j++;
        }

        if (i == n && j == t.size())
            return 0;
        if (i == n)
            return -1;
        if (j == t.size())
            return 1;
        if (s[i] < t[j])return -1;
        return 1;
    }
    bool is_prefix(int ind, const string& t) {
        bool ok = true;
        int i = suf[ind], j = 0;
        while (i < n && j < t.size() && ok) {
            ok &= s[i] == t[j];
            i++, j++;
        }
        return ok;
    }
    int find_string(const string& t) {
        int l = 0, r = n - 1, ans = -1;
        while (l <= r) {
            int mid = (l + r) / 2;
            if (comp(suf[mid], t) >= 0) {
                ans = mid;
                r = mid - 1;
            }
            else {
                l = mid + 1;
            }
        }
        if (ans == -1)
            return -1;
        return is_prefix(ans, t) ? ans : -1;
    }
    int count_string(const string& t) {
        int index = find_string(t);
        if (index == -1)
            return 0;

        int l = index, r = n - 1, ans = -1;
        while (l <= r) {
            int mid = (l + r) / 2;
            if (is_prefix(mid, t)) {
                ans = mid;
                l = mid + 1;
            }
            else {
                r = mid - 1;
            }
        }
        return ans - index + 1;
    }
    void calc_lcp() {
        lcp.resize(n);
        int k = 0;
        for (int i = 0; i < n - 1; i++) {
            int pi = classes[i];
            int j = suf[pi - 1];
            while (s[i + k] == s[j + k])k++;
            lcp[pi] = k;
            k = max(k - 1, 0);
        }
    }
};
long long calc_distinct_substring(string& s) {
    SuffixArray sa(s);
    long long ans = 0;
    for (int i = 1; i < sa.lcp.size(); i++) {
        int sz = sa.n - sa.suf[i] - 1;
        ans += (sz - sa.lcp[i]);
    }
    return ans;
}
pair<int, int> lcs(string& s, string& t) {
    int n = s.size(), m = t.size();
    string y = s + '$' + t;
    SuffixArray sa(y); // 0..n-1 -> s, n+1..n+m -> t

    // for (int i = 0; i < sa.n; i++) {
    //     cout << sa.suf[i] << " " << sa.lcp[i] << " " << sa.s.substr(sa.suf[i]) << "\n";
    // }

    int maximum = 0, index = 0;
    for (int i = 0; i + 1 <= n + m + 1; i++) {
        if (sa.suf[i] <= n - 1 && n + 1 <= sa.suf[i + 1]) {
            if (maximum < sa.lcp[i + 1]) {
                index = sa.suf[i];
                maximum = sa.lcp[i + 1];
            }
        }
        else if (sa.suf[i + 1] <= n - 1 && n + 1 <= sa.suf[i]) {
            if (maximum < sa.lcp[i + 1]) {
                index = sa.suf[i + 1];
                maximum = sa.lcp[i + 1];
            }
        }
    }
    return {index, maximum};
}

// ==========================================================================================

struct aho_corasick {
    const int K = 26;
    struct vertex {
        char ch;
        bool leaf = false;
        int next[K], parent, fail, cnt;

        vertex(int _p = -1, char _ch = '$') : parent(_p), ch(_ch) {
            fail = cnt = 0;
            memset(next, -1, sizeof next);
        }
    };

    int time = -1;
    vector<vertex> tree;
    vector<int> in, out;
    aho_corasick() { tree.emplace_back(); }

    void insert(const string &s, int idx) {
        int u = 0;
        for (auto &c: s) {
            int v = c - 'a';
            if (tree[u].next[v] == -1) {
                tree[u].next[v] = tree.size();
                tree.emplace_back(u, c);
            }
            u = tree[u].next[v];
        }
        tree[u].cnt++;
        tree[u].leaf = true;
    }
    void dfs(int u) {
        in[u] = ++time;
        for (int i = 0; i < K; i++) {
            if (tree[u].next[i] != -1) {
                dfs(tree[u].next[i]);
            }
        }
        out[u] = time;
    }
    bool is_ancestor(int u, int v) {
        return in[u] <= in[v] && out[v] <= out[u];
    }
    void build_aho() {
        in.resize(tree.size());
        out.resize(tree.size());
        dfs(0);

        queue<int> q;
        q.push(0);

        auto compute_failure = [&](int u) -> int {
            if (u == 0 || tree[u].parent == 0)
                return 0;

            int x = tree[u].ch - 'a';
            while (u) {
                int p = tree[u].parent;
                int f = tree[p].fail;
                if (tree[f].next[x] != -1) {
                    return tree[f].next[x];
                }
                u = f;
            }
            if (tree[u].next[x] != -1)
                return tree[u].next[x];
            return u;
        };

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            tree[u].fail = compute_failure(u);
            if (~tree[u].parent)
                tree[u].cnt += tree[ tree[u].parent ].cnt;
            for (int i = 0; i < K; i++) {
                if (tree[u].next[i] != -1) {
                    q.push(tree[u].next[i]);
                }
            }
        }
    }
    void print() {
        // tree
        for (int i = 0; i < tree.size(); i++) {
            for (int j = 0; j < K; j++) {
                if (tree[i].next[j] != -1) {
                    cout << i << " ";
                    cout << tree[i].next[j] << " ";
                    cout << char(j + 'a') << "\n";

                }
            }
        }
        // failure
        for (int i = 0; i < tree.size(); i++) {
            cout << i << " ";
            cout << tree[i].fail << "\n";
        }
    }

    int fail(int u) { return tree[u].fail; }    
    int next(int u, char c) { return tree[u].next[c]; }
    vertex operator[](int u) { return tree[u]; }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
