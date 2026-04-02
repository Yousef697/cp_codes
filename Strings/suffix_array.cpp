#include <bits/stdc++.h>

using namespace std;
using ll = long long;

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

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    return 0;
}
