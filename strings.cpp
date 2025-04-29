#include <bits/stdc++.h>

using namespace std;

// KMP, String Matching
/*
    prefix: any sub string that starts from index 0
    suffix: any sub string that ends from index size-1
    proper prefix: any prefix that its size is not equal to the main string size

    compute_prefixes function (failure function) compute the proper prefixes of a string such that
    ret[i] = max(proper prefix: prroper prefix = suffix )
*/
vector<int> compute_prefixes(vector<int> pat)
{

    int m = pat.size();

    int i;
    int k;

    vector<int> prefixes(m);

    for (i = 1, k = 0; i < m; i++)
    {

        while (k > 0 && pat[k] != pat[i])
        {
            k = prefixes[k - 1]; // return to the promising position
        }

        if (pat[i] == pat[k])
            prefixes[i] = ++k;
        else
            prefixes[i] = k;
    }

    return prefixes;
}
int KMP(vector<int> s, vector<int> pat)
{

    int n = s.size(), m = pat.size(), ans = 0;

    vector<int> prefixes = compute_prefixes(pat);

    int i; // for s
    int k; // for pat

    for (i = 0, k = 0; i < n; i++)
    {

        while (k > 0 && pat[k] != s[i])
        {
            k = prefixes[k - 1]; // back to the promising position
        }

        if (s[i] == pat[k])
            k++;

        if (k == m)
            ans++, k = prefixes[k - 1]; // make k fail
    }

    return ans;
}

// ==========================================================================================

// Z-Algorithm
vector<int> z_algorithm(const string& s)
{
    int n = s.size();
    vector<int> z(n);
    int l = 0, r = 0;

    for (int i = 1; i < n; i++)
    {
        if (i < r)
            z[i] = min(r - i, z[i - l]);

        while (i + z[i] < n && s[i + z[i]] == s[z[i]])
            z[i]++;

        if (i + z[i] > r)
            l = i, r = i + z[i];
    }

    return z;
}
vector<int> find_pattern(const string& s, const string& t)
{
    string A = s + '#' + t;
    int n = s.size(), m = t.size();
    vector<int> Z = z_algorithm(A);

    vector<int> ans;
    for (int i = 0; i < A.size(); i++)
    {
        if (Z[i + n + 1] == n)
            ans.push_back(i);
    }

    return ans;
}
int distinct_substrings(const string& s)
{
    int n = s.size(), m = 0, ans = 0;
    string t = "";

    for (const char& c : s)
    {
        m++;
        t += c;
        int mx = 0;
        reverse(t.begin(), t.end());

        vector<int> Z = z_algorithm(t);
        for (const int i : Z)
            mx = max(mx, i);
        ans += m - mx;

        reverse(t.begin(), t.end());
    }

    return ans;
}
string period(const string& s)
{
    int n = s.size(), ans = 0;
    vector<int> Z = z_algorithm(s);

    for (int i = 1; i <= n; i++)
    {
        if (n % i)
            continue;

        if (i + Z[i] == n)
        {
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

struct Hash
{
    int mod, base, n;
    vector<int> prf, suf, pows, invs;
    vector<int> primesMod = {1000000007, 1000000009, 1000000021, 1000000033, 1000000087, 1000000093, 1000000097, 1000000103, 1000000123, 1000000181, 1000000207, 1000000223, 1000000241, 1000000271, 1000000289, 1000000297, 1000000321, 1000000349, 1000000363, 1000000403, 1000000409, 1000000411, 1000000427, 1000000433, 1000000439, 1000000447, 1000000453, 1000000459, 1000000483, 1000000513, 1000000531, 1000000579, 1000000607, 1000000613, 1000000637, 1000000663, 1000000711, 1000000753, 1000000787, 1000000801, 1000000829, 1000000861, 1000000871, 1000000891, 1000000901, 1000000919, 1000000931, 1000000933, 1000000993, 1000001011};

    int fpow(int x, int p)
    {
        int ans = 1;
        while (p)
        {
            if (p & 1)
                ans = (1LL * ans * x) % mod;
            p /= 2;
            x = (1LL * x * x) % mod;
        }
        return ans;
    }

    Hash(const string& s)
    {
        int r = rand(0, 49);
        mod = primesMod[r];
        base = rand(2, 199);

        n = s.size();
        pows.assign(n + 5, 1);
        invs.assign(n + 5, 1);

        for (int i = 1; i < n + 5; i++)
        {
            pows[i] = (1LL * pows[i - 1] * base) % mod;
            invs[i] = fpow(pows[i], mod - 2);
        }

        int value = 0, p = 0;
        prf.push_back(0);
        for (int i = 0; i < n; i++)
        {
            char c = s[i];
            value = (value + 1LL * (c - 'a' + 1) * pows[p++]) % mod;
            prf.push_back(value);
        }
        prf.push_back(0);

        value = 0, p = 0;
        suf.push_back(0);
        for (int i = n - 1; i >= 0; i--)
        {
            char c = s[i];
            value = (value + 1LL * (c - 'a' + 1) * pows[p++]) % mod;
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
        }
        else /// suffix
        {
            int ll = l, rr = r;
            l = n - rr + 1, r = n - ll + 1;
            cand = (suf[r] - suf[l - 1]) % mod;
            cand = (cand + mod) % mod;
            cand = (1LL * cand * invs[l]) % mod;
        }
        return cand;
    }

    bool palindrome(int l, int r)
    {
        return substring(l, r, 0) == substring(l, r, 1);
    }
};

// ==========================================================================================

// Trie
const int K = 26;
struct Trie
{
    vector<int> leaf;
    vector<vector<int>> tree;

    Trie()
    {
        leaf.push_back(0);
        tree.push_back(vector<int>(26, 0));
    }

    void add(const string& s)
    {
        int v = 0;
        for (const char& c : s)
        {
            if (tree[v][c - 'a'] == 0)
            {
                tree[v][c - 'a'] = tree.size();
                tree.push_back(vector<int>(26, 0));
                leaf.push_back(0);
            }
            v = tree[v][c - 'a'];
        }
        leaf[v] = true;
    }

    int serach_word(const string& s)
    {
        int v = 0;
        for (const char& c : s)
        {
            if (tree[v][c - 'a'] == 0)
                return -1;

            v = tree[v][c - 'a'];
        }
        return (leaf[v] == 1 ? 1 : -1);
    }

    int serach_prefix(const string& s)
    {
        int v = 0;
        for (const char& c : s)
        {
            if (tree[v][c - 'a'] == 0)
                return -1;

            v = tree[v][c - 'a'];
        }
        return v;
    }
};

// ==========================================================================================

// Recursive Trie
const int K = 2;
struct Trie
{
    vector<int> leaf;
    vector<int> cnts;
    vector<vector<int>> tree;

    Trie()
    {
        leaf.push_back(0);
        cnts.push_back(0);
        tree.push_back(vector<int>(2, 0));
    }

    void insert(const string& s, int v = 0, int i = 0)
    {
        if (i == s.size())
        {
            cnts[v]++;
            leaf[v] = 1;
            return;
        }

        if (!tree[v][s[i] - '0'])
        {
            tree[v][s[i] - '0'] = tree.size();

            leaf.push_back(0);
            cnts.push_back(0);
            tree.push_back(vector<int>(2, 0));
        }

        int next = tree[v][s[i] - '0'];
        insert(s, next, i + 1);

        cnts[v]++;
    }

    void erase(const string& s, int v = 0, int i = 0)
    {
        if (i == s.size())
        {
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

    int serach_word(const string& s, int v = 0, int i = 0)
    {
        if (i == s.size())
            return leaf[v];

        if (!tree[v][s[i] - '0'])
            return -1;
        return serach_word(s, tree[v][s[i] - '0'], i + 1);
    }

    int serach_prefix(const string& s, int v = 0, int i = 0)
    {
        if (i == s.size())
            return 1;

        if (!tree[v][s[i] - '0'])
            return -1;
        return serach_prefix(s, tree[v][s[i] - '0'], i + 1);
    }
};

// ==========================================================================================

// Manacher
vector<int> manacher(const string& s)
{
    string t = "";
    for (const char& c : s)
        t += '#', t += c;
    t += '#';
    int n = t.size();
    t = '!' + t + '*';

    int l = 0, r = 1;
    vector<int> p(n + 2);
    for (int i = 1; i <= n; i++)
    {
        p[i] = min(r - i, p[ l + (r - i) ]);
        p[i] = max(0, p[i]);

        while (t[i - p[i]] == t[i + p[i]])
            p[i]++;

        if (i + p[i] > r)
            l = i - p[i], r = i + p[i];
    }

    for (int &i : p)
        i--;

    return vector<int>(p.begin() + 1, p.end() - 1);
}

// ==========================================================================================

int32_t main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}