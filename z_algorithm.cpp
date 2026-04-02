#include <bits/stdc++.h>

using namespace std;
using ll = long long;

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

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    return 0;
}
