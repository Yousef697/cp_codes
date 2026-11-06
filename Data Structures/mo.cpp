#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using ld = long double;

struct MO {
    struct Query {
        int l, r, id, b;
        bool operator <(const Query &other) const {
            if (b != other.b) return b < other.b;
            return b & 1 ? r > other.r : r < other.r;
        }
    };

    int n, sq;
    vector<int> values;
    vector<Query> queries;

    MO(const vector<int> &values): n(values.size()), sq(sqrt(values.size()) + 1), values(values) {}

    void add_query(int l, int r, int id) {
        queries.push_back({l, r, id, l / sq});
    }

    vector<int> solve_queries() {
        sort(queries.begin(), queries.end());

        int l = 0, r = -1, cnt = 0;
        vector<int> ans(queries.size());
        vector<int> freq(*max_element(values.begin(), values.end()) + 1);

        auto push = [&](int idx) -> void {
            int val = values[idx];
            freq[val]++;
            if (freq[val] == 1)
                cnt++;
        };
        auto pop = [&](int idx) -> void {
            int val = values[idx];
            freq[val]--;
            if (freq[val] == 0)
                cnt--;
        };
        auto solve = [&]() -> int {
            return cnt;
        };
        for (auto &q: queries) {
            while (r < q.r) push(++r);
            while (l > q.l) push(--l);
            while (l < q.l) pop(l++);
            while (r > q.r) pop(r--);
            ans[q.id] = solve();
        }
        return ans;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;

    vector<int> a(n);
    for (auto& i : a)
        cin >> i;
    MO mo(a);

    int q, c = 0;
    cin >> q;
    while (q--) {
        int l, r;
        cin >> l >> r;
        l--, r--;
        mo.add_query(l, r, c);
        c++;
    }
    auto ans = mo.solve_queries();
    for (auto& i : ans)
        cout << i << "\n";

    return 0;
}
