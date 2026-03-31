#include <bits/stdc++.h>

using namespace std;
using ll = long long;

// Ordered Set / Ordered Multiset
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
template<class T> using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template<class T> using ordered_multiset = tree<T, null_type, less_equal<T>, rb_tree_tag, tree_order_statistics_node_update>;

// Compress Numbers
int compress_vector(vector<int>& ranges) {
    vector<int> c;
    for (auto& i : ranges) {
        for (int k = -1; k <= 1; k++) {
            c.emplace_back(i + k);
        }
    }
    sort(c.begin(), c.end());
    c.erase(unique(c.begin(), c.end()), c.end());
    for (auto& i : ranges) {
        i = lower_bound(c.begin(), c.end(), i) - c.begin();
    }
    return c.back(); // max element if needed
}

// Compress Ranges
int compress_ranges(vector<pair<int, int>>& ranges) {
    vector<int> c;
    for (auto& [i, j] : ranges) {
        for (int k = -1; k <= 1; k++) {
            c.emplace_back(i + k);
            c.emplace_back(j + k);
        }
    }
    sort(c.begin(), c.end());
    c.erase(unique(c.begin(), c.end()), c.end());
    for (auto& [i, j] : ranges) {
        i = lower_bound(c.begin(), c.end(), i) - c.begin();
        j = lower_bound(c.begin(), c.end(), j) - c.begin();
    }
    return c.back(); // max element if needed
}

// Union Ranges
void ranges_union(vector<pair<int, int>>& ranges) {

    vector<pair<int, int>> cur;
    sort(ranges.begin(), ranges.end());

    cur.emplace_back(ranges[0]);
    for (int i = 1; i < ranges.size(); i++) {
        if (cur.back().first == ranges[i].first)
            cur.back().second = ranges[i].second;
        else if (cur.back().second < ranges[i].first)
            cur.push_back(ranges[i]);
        else
            cur.back().second = max(cur.back().second, ranges[i].second);
    }
    
    ranges = cur;
}


int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
