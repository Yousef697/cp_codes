#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using ld = long double;

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
    return *max_element(ranges.begin(), ranges.end()); // max element if needed
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
        if (cur.back().second < ranges[i].first)
            cur.push_back(ranges[i]);
        else
            cur.back().second = max(cur.back().second, ranges[i].second);
    }
    
    ranges = cur;
}

// fraction struct
struct fraction {
    ll a, b;
 
    fraction() : a(0), b(1) {}
    fraction(ll x) : a(x), b(1) {}
    fraction(ll x, ll y) : a(x), b(y) { norm(); }
 
    ll gcd(ll x, ll y) { return y == 0 ? x : gcd(y, x % y); }
    void norm() { ll g = gcd(a, b); a /= g, b /= g; }
};
fraction add(const fraction& n, const fraction& m) { return {n.a * m.b + n.b * m.a, n.b * m.b}; }
fraction sub(const fraction& n, const fraction& m) { return {n.a * m.b - n.b * m.a, n.b * m.b}; }
fraction mul(const fraction& n, const fraction& m) { return {n.a * m.a, n.b * m.b}; }
fraction div(const fraction& n, const fraction& m) { return {n.a * m.b,  n.b * m.a}; }
fraction abs(const fraction& n) { return {abs(n.a), abs(n.b)}; }
int cmp(const fraction& n, const fraction& m) {
    if (n.a == m.a && n.b == m.b) return 0;
    return n.a * m.b < m.a * n.b ? -1 : 1;
}
void println(const fraction& n) { cout << n.a << "/" << n.b << "\n"; }
void print(fraction& n) { cout << n.a << "/" << n.b; }

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
