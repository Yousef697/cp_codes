#include <bits/stdc++.h>

using namespace std;

// Segment Tree
// For Min
template <typename T>
struct SegmentTree
{
    int size;
    vector<T> v, seg;

    SegmentTree() {}

    SegmentTree(vector<T> vec)
    {
        v = vec;
        size = 1;
        int n = vec.size() - 1;
        while (size < n)
            size *= 2;
        seg.assign(2 * size + 2, 0);
        build();
    }

    T merge(T a, T b)
    {
        return min(a, b);
    }

    void build(int n, int l, int r)
    {
        if (l == r)
        {
            if (l < v.size())
                seg[n] = v[l];
            return;
        }

        build(2 * n, l, (l + r) / 2);
        build(2 * n + 1, (l + r) / 2 + 1, r);

        seg[n] = merge(seg[2 * n], seg[2 * n + 1]);
    }

    void update(int n, int l, int r, int ind, int val)
    {
        if (ind < l || ind > r)
            return;

        if (l == r)
        {
            if (l < v.size())
                seg[n] = val, v[l] = val;
            return;
        }

        update(2 * n, l, (l + r) / 2, ind, val);
        update(2 * n + 1, (l + r) / 2 + 1, r, ind, val);

        seg[n] = merge(seg[2 * n], seg[2 * n + 1]);
    }

    T query(int n, int l, int r, int lq, int rq)
    {
        if (r < lq || l > rq)
            return INT_MAX;

        if (lq <= l && r <= rq)
            return seg[n];

        T a = query(2 * n, l, (l + r) / 2, lq, rq);
        T b = query(2 * n + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge(a, b);
    }

    void build() { build(1, 1, size); }

    void update(int ind, int val) { update(1, 1, size, ind, val); }

    T query(int l, int r) { return query(1, 1, size, l, r); }
};

// For Max
template <typename T>
struct SegmentTree
{
    int size;
    vector<T> v, seg;

    SegmentTree() {}

    SegmentTree(vector<T> vec)
    {
        v = vec;
        size = 1;
        int n = vec.size() - 1;
        while (size < n)
            size *= 2;
        seg.assign(2 * size + 2, 0);
        build();
    }

    T merge(T a, T b)
    {
        return max(a, b);
    }

    void build(int n, int l, int r)
    {
        if (l == r)
        {
            if (l < v.size())
                seg[n] = v[l];
            return;
        }

        build(2 * n, l, (l + r) / 2);
        build(2 * n + 1, (l + r) / 2 + 1, r);

        seg[n] = merge(seg[2 * n], seg[2 * n + 1]);
    }

    void update(int n, int l, int r, int ind, int val)
    {
        if (ind < l || ind > r)
            return;

        if (l == r)
        {
            if (l < v.size())
                seg[n] = val, v[l] = val;
            return;
        }

        update(2 * n, l, (l + r) / 2, ind, val);
        update(2 * n + 1, (l + r) / 2 + 1, r, ind, val);

        seg[n] = merge(seg[2 * n], seg[2 * n + 1]);
    }

    T query(int n, int l, int r, int lq, int rq)
    {
        if (r < lq || l > rq)
            return -1;

        if (lq <= l && r <= rq)
            return seg[n];

        T a = query(2 * n, l, (l + r) / 2, lq, rq);
        T b = query(2 * n + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge(a, b);
    }

    void build() { build(1, 1, size); }

    void update(int ind, int val) { update(1, 1, size, ind, val); }

    T query(int l, int r) { return query(1, 1, size, l, r); }
};

// For Sum
template <typename T>
struct SegmentTree
{
    int size;
    vector<T> v, seg;

    SegmentTree() {}

    SegmentTree(vector<T> vec)
    {
        v = vec;
        size = 1;
        int n = vec.size() - 1;
        while (size < n)
            size *= 2;
        seg.assign(2 * size + 2, 0);
        build();
    }

    T merge(T a, T b)
    {
        return a + b;
    }

    void build(int n, int l, int r)
    {
        if (l == r)
        {
            if (l < v.size())
                seg[n] = v[l];
            return;
        }

        build(2 * n, l, (l + r) / 2);
        build(2 * n + 1, (l + r) / 2 + 1, r);

        seg[n] = merge(seg[2 * n], seg[2 * n + 1]);
    }

    void update(int n, int l, int r, int ind, int val)
    {
        if (ind < l || ind > r)
            return;

        if (l == r)
        {
            if (l < v.size())
                seg[n] = val, v[l] = val;
            return;
        }

        update(2 * n, l, (l + r) / 2, ind, val);
        update(2 * n + 1, (l + r) / 2 + 1, r, ind, val);

        seg[n] = merge(seg[2 * n], seg[2 * n + 1]);
    }

    T query(int n, int l, int r, int lq, int rq)
    {
        if (r < lq || l > rq)
            return 0;

        if (lq <= l && r <= rq)
            return seg[n];

        T a = query(2 * n, l, (l + r) / 2, lq, rq);
        T b = query(2 * n + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge(a, b);
    }

    void build() { build(1, 1, size); }

    void update(int ind, int val) { update(1, 1, size, ind, val); }

    T query(int l, int r) { return query(1, 1, size, l, r); }
};

// MinMax Segment Tree
template <typename T>
struct SegmentTree
{
    int size;
    vector<T> v, seg_min, seg_max;

    SegmentTree() {}

    SegmentTree(vector<T> vec)
    {
        v = vec;
        size = 1;
        int n = vec.size() - 1;
        while (size < n)
            size *= 2;
        seg_min.assign(2 * size + 2, 2e9);
        seg_max.assign(2 * size + 2, -2e9);

        build();
    }

    T merge_min(T a, T b) { return min(a, b); }
    T merge_max(T a, T b) { return max(a, b); }

    void build(int n, int l, int r)
    {
        if (l == r)
        {
            if (l < v.size())
                seg_min[n] = seg_max[n] = v[l];
            return;
        }

        build(2 * n, l, (l + r) / 2);
        build(2 * n + 1, (l + r) / 2 + 1, r);

        seg_min[n] = merge_min(seg_min[2 * n], seg_min[2 * n + 1]);
        seg_max[n] = merge_max(seg_max[2 * n], seg_max[2 * n + 1]);
    }

    void update(int n, int l, int r, int ind, int val)
    {
        if (ind < l || ind > r)
            return;

        if (l == r)
        {
            if (l < v.size())
                seg_min[n] = seg_max[n] = v[l];
            return;
        }

        update(2 * n, l, (l + r) / 2, ind, val);
        update(2 * n + 1, (l + r) / 2 + 1, r, ind, val);

        seg_min[n] = merge_min(seg_min[2 * n], seg_min[2 * n + 1]);
        seg_max[n] = merge_max(seg_max[2 * n], seg_max[2 * n + 1]);
    }

    T query_min(int n, int l, int r, int lq, int rq)
    {
        if (r < lq || l > rq)
            return 2e9;

        if (lq <= l && r <= rq)
            return seg_min[n];

        T a = query_min(2 * n, l, (l + r) / 2, lq, rq);
        T b = query_min(2 * n + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge_min(a, b);
    }

    T query_max(int n, int l, int r, int lq, int rq)
    {
        if (r < lq || l > rq)
            return -2e9;

        if (lq <= l && r <= rq)
            return seg_max[n];

        T a = query_max(2 * n, l, (l + r) / 2, lq, rq);
        T b = query_max(2 * n + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge_max(a, b);
    }

    void build() { build(1, 1, size); }

    void update(int ind, int val) { update(1, 1, size, ind, val); }

    T query_min(int l, int r) { return query_min(1, 1, size, l, r); }

    T query_max(int l, int r) { return query_max(1, 1, size, l, r); }
};

// ==========================================================================================

// Lazy Propagation, Lazy Segment Tree
// Lazy For Min
struct LazySegmentTree
{
    int size;
    vector<long long> v, seg, lazy;

    LazySegmentTree(vector<long long> vec)
    {
        v = vec;
        size = 1;
        int n = vec.size() - 1;
        while (size < n)
            size *= 2;
        seg.assign(2 * size + 2, 0);
        lazy.assign(2 * size + 2, 0);
    }

    void process(int n, int l, int r)
    {
        if (lazy[n] == 0)
            return;

        seg[n] += lazy[n];
        if (l != r)
            lazy[2 * n] += lazy[n], lazy[2 * n + 1] += lazy[n];
        lazy[n] = 0;
    }

    void build(int n, int l, int r)
    {
        process(n, l, r);

        if (l == r)
        {
            if (l < v.size())
                seg[n] = v[l];
            return;
        }

        build(2 * n, l, (l + r) / 2);
        build(2 * n + 1, (l + r) / 2 + 1, r);

        seg[n] = min(seg[2 * n], seg[2 * n + 1]);
    }

    void update_range(int n, int l, int r, int lq, int rq, int x)
    {
        process(n, l, r);

        if (r < lq || rq < l)
            return;

        if (lq <= l && r <= rq)
        {
            seg[n] += x;
            if (l != r)
                lazy[2 * n] += x, lazy[2 * n + 1] += x;
            return;
        }

        update_range(2 * n, l, (l + r) / 2, lq, rq, x);
        update_range(2 * n + 1, (l + r) / 2 + 1, r, lq, rq, x);

        seg[n] = min(seg[2 * n], seg[2 * n + 1]);
    }

    long long query(int n, int l, int r, int lq, int rq)
    {
        process(n, l, r);

        if (r < lq || l > rq)
            return LLONG_MAX;

        if (lq <= l && r <= rq)
            return seg[n];

        int a = query(2 * n, l, (l + r) / 2, lq, rq);
        int b = query(2 * n + 1, (l + r) / 2 + 1, r, lq, rq);

        return min(a, b);
    }

    void build() { build(1, 1, size); }

    void update_range(int l, int r, int x) { update_range(1, 1, size, l, r, x); }

    long long query(int l, int r) { return query(1, 1, size, l, r); }
};

// Lazy For Sum
struct LazySegmentTree
{
    int size;
    vector<long long> v, seg, lazy;

    LazySegmentTree(vector<long long> vec)
    {
        v = vec;
        size = 1;
        int n = vec.size() - 1;
        while (size < n)
            size *= 2;
        seg.assign(2 * size + 2, 0);
        lazy.assign(2 * size + 2, 0);
    }

    void process(int n, int l, int r)
    {
        if (lazy[n] == 0)
            return;

        seg[n] += lazy[n] * (r - l + 1);
        if (l != r)
            lazy[2 * n] += lazy[n], lazy[2 * n + 1] += lazy[n];
        lazy[n] = 0;
    }

    void build(int n, int l, int r)
    {
        process(n, l, r);

        if (l == r)
        {
            if (l < v.size())
                seg[n] = v[l];
            return;
        }

        build(2 * n, l, (l + r) / 2);
        build(2 * n + 1, (l + r) / 2 + 1, r);

        seg[n] = seg[2 * n] + seg[2 * n + 1];
    }

    void update_range(int n, int l, int r, int lq, int rq, long long x)
    {
        process(n, l, r);

        if (r < lq || rq < l)
            return;

        if (lq <= l && r <= rq)
        {
            seg[n] += x * (r - l + 1);
            if (l != r)
                lazy[2 * n] += x, lazy[2 * n + 1] += x;
            return;
        }

        update_range(2 * n, l, (l + r) / 2, lq, rq, x);
        update_range(2 * n + 1, (l + r) / 2 + 1, r, lq, rq, x);

        seg[n] = seg[2 * n] + seg[2 * n + 1];
    }

    long long query(int n, int l, int r, int lq, int rq)
    {
        process(n, l, r);

        if (r < lq || l > rq)
            return 0;

        if (lq <= l && r <= rq)
            return seg[n];

        long long a = query(2 * n, l, (l + r) / 2, lq, rq);
        long long b = query(2 * n + 1, (l + r) / 2 + 1, r, lq, rq);

        return a + b;
    }

    void build() { build(1, 1, size); }

    void update_range(int l, int r, long long x) { update_range(1, 1, size, l, r, x); }

    long long query(int l, int r) { return query(1, 1, size, l, r); }
};

// ==========================================================================================

// Maximum Sum Range Query
struct node
{
    long long prefix, middle, suffix, sum;

    node(long long a, long long b, long long c, long long d)
    {
        prefix = a, middle = b, suffix = c, sum = d;
    }
};
struct MaximumSum
{
    int size = 1;
    vector<node> seg;

    MaximumSum(int n)
    {
        while (size < n)
            size *= 2;
        seg.assign(2 * size + 4, node(0, 0, 0, 0));
    }

    node merge(int n, node a, node b)
    {
        if (a.prefix == -OO)
            return b;
        if (b.prefix == -OO)
            return a;

        node ans(0, 0, 0, 0);

        ans.prefix = max({a.prefix, a.sum + b.prefix, a.sum + b.sum});

        ans.middle = max({a.middle, b.middle, a.suffix + b.prefix, a.sum + b.sum});

        ans.suffix = max({b.suffix, b.sum + a.suffix, a.sum + b.sum});

        ans.sum = a.sum + b.sum;

        return ans;
    }
    void build(vector<long long> &v, int n, int l, int r)
    {
        if (l == r)
        {
            if (l < v.size())
                seg[n] = node(v[l], v[l], v[l], v[l]);
            return;
        }

        build(v, 2 * n, l, (l + r) / 2);
        build(v, 2 * n + 1, (l + r) / 2 + 1, r);

        seg[n] = merge(n, seg[2 * n], seg[2 * n + 1]);
    }
    void update(int n, int l, int r, int ind, long long val)
    {
        if (ind < l || ind > r)
            return;

        if (l == r)
        {
            seg[n].prefix = val;
            seg[n].middle = val;
            seg[n].suffix = val;
            seg[n].sum = val;
            return;
        }

        update(2 * n, l, (l + r) / 2, ind, val);
        update(2 * n + 1, (l + r) / 2 + 1, r, ind, val);

        seg[n] = merge(n, seg[2 * n], seg[2 * n + 1]);
    }
    node query(int n, int l, int r, int lq, int rq)
    {
        if (l > rq || lq > r)
            return node(-OO, -OO, -OO, -OO);

        if (lq <= l && r <= rq)
            return seg[n];

        node a = query(2 * n, l, (l + r) / 2, lq, rq);
        node b = query(2 * n + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge(n, a, b);
    }

    void build(vector<long long> &v) { build(v, 1, 1, size); }

    void update(int i, long long val) { update(1, 1, size, i, val); }

    long long query(int l, int r)
    {
        node ans = query(1, 1, size, l, r);
        return max({ans.prefix, ans.middle, ans.suffix});
    }
};

// ==========================================================================================

// Maximum Alternating Subsequence Sum
struct node
{
    long long max_pp; // max of start with + and end with +
    long long max_pn; // max of start with + and end with -

    long long min_pp; // min of start with + and end with +
    long long min_pn; // min of start with + and end with -

    node(long long a, long long b)
    {
        max_pp = a;
        min_pp = a;

        max_pn = b;
        min_pn = b;
    }
};
struct MaximumAlternatingSubsequenceSum
{
    int size = 1;
    vector<long long> v;
    vector<node> seg;

    MaximumAlternatingSubsequenceSum(vector<long long> &vec)
    {
        v = vec;
        int n = vec.size() - 1;
        while (size < n)
            size *= 2;
        seg.assign(2 * size + 4, node(0, 0));
        build();
    }

    node merge(node a, node b)
    {
        node ans(0, 0);

        ans.max_pp = max({a.max_pp, b.max_pp, a.max_pp - b.min_pn, a.max_pn + b.max_pp});
        ans.max_pn = max({a.max_pn, b.max_pn, a.max_pn + b.max_pn, a.max_pp - b.min_pp});

        ans.min_pp = min({a.min_pp, b.min_pp, a.min_pp - b.max_pn, a.min_pn + b.min_pp});
        ans.min_pn = min({a.min_pn, b.min_pn, a.min_pn + b.min_pn, a.min_pp - b.max_pp});

        return ans;
    }

    void build(int n, int l, int r)
    {
        if (l == r)
        {
            if (l < v.size())
                seg[n] = node(v[l], 0);
            return;
        }

        build(2 * n, l, (l + r) / 2);
        build(2 * n + 1, (l + r) / 2 + 1, r);

        seg[n] = merge(seg[2 * n], seg[2 * n + 1]);
    }

    void update(int n, int l, int r, int ind, long long val)
    {
        if (ind < l || ind > r)
            return;

        if (l == r)
        {
            v[ind] = val;
            seg[n] = node(val, 0);

            return;
        }

        update(2 * n, l, (l + r) / 2, ind, val);
        update(2 * n + 1, (l + r) / 2 + 1, r, ind, val);

        seg[n] = merge(seg[2 * n], seg[2 * n + 1]);
    }

    node query(int n, int l, int r, int lq, int rq)
    {
        if (l > rq || lq > r)
            return node(0, 0);

        if (lq <= l && r <= rq)
            return seg[n];

        node a = query(2 * n, l, (l + r) / 2, lq, rq);
        node b = query(2 * n + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge(a, b);
    }

    void build() { build(1, 1, size); }

    void update(int i, long long val) { update(1, 1, size, i, val); }

    long long query(int l, int r)
    {
        node ans = query(1, 1, size, l, r);
        return max({ans.max_pp, ans.max_pn});
    }
};

// ==========================================================================================

// Merge Sort Tree
struct MergeSortTree
{
    int size;
    vector<int> v;
    vector<vector<int>> seg;

    MergeSortTree() {}

    MergeSortTree(vector<int> &vec)
    {
        v = vec;
        size = 1;
        int n = vec.size() - 1;
        while (size < n)
            size *= 2;
        seg.assign(2 * size + 2, vector<int>{});

        build();
    }

    vector<int> merge(vector<int> &a, vector<int> &b)
    {
        vector<int> ans;
        int x = a.size(), y = b.size(), i = 0, j = 0;

        while (i < x && j < y)
        {
            if (a[i] < b[j])
                ans.emplace_back(a[i++]);
            else
                ans.emplace_back(b[j++]);
        }

        while (i < x)
            ans.emplace_back(a[i++]);

        while (j < y)
            ans.emplace_back(b[j++]);

        return ans;
    }

    void build(int n, int l, int r)
    {
        if (l == r)
        {
            if (l < v.size())
                seg[n] = vector<int>{v[l]};
            return;
        }

        build(2 * n, l, (l + r) / 2);
        build(2 * n + 1, (l + r) / 2 + 1, r);

        seg[n] = merge(seg[2 * n], seg[2 * n + 1]);
    }

    int greater(int n, int l, int r, int lq, int rq, int val)
    {
        if (r < lq || l > rq)
            return 0;

        if (lq <= l && r <= rq)
        {
            int ind = upper_bound(seg[n].begin(), seg[n].end(), val) - seg[n].begin();
            return seg[n].size() - ind;
        }

        int a = greater(2 * n, l, (l + r) / 2, lq, rq, val);
        int b = greater(2 * n + 1, (l + r) / 2 + 1, r, lq, rq, val);

        return (a + b);
    }

    int less(int n, int l, int r, int lq, int rq, int val)
    {
        if (r < lq || l > rq)
            return 0;

        if (lq <= l && r <= rq)
        {
            int ind = lower_bound(seg[n].begin(), seg[n].end(), val) - seg[n].begin() - 1;
            return ind + 1;
        }

        int a = less(2 * n, l, (l + r) / 2, lq, rq, val);
        int b = less(2 * n + 1, (l + r) / 2 + 1, r, lq, rq, val);

        return (a + b);
    }

    void build() { build(1, 1, size); }

    int less(int l, int r, int val) { return less(1, 1, size, l, r, val); }

    int greater(int l, int r, int val) { return greater(1, 1, size, l, r, val); }
};

int32_t main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}