#include <bits/stdc++.h>
#define int long long

using namespace std;

// Iterative Segment Tree
struct Segment_Tree {
    using T = int;

    int size, neutral;
    vector<int> seg;

    int merge(int a, int b) { return max(a, b); }

    Segment_Tree() {
    }

    Segment_Tree(int n) {
        size = 1, neutral = 0;
        while (size < n) size <<= 1;
        seg.assign(2 * size + 2, 0);
    }

    Segment_Tree(const vector<int> &vec) {
        int n = vec.size();
        size = 1, neutral = 0;
        while (size < n) size <<= 1;
        seg.assign(2 * size + 2, 0);
        for (int i = 0; i < n; ++i)
            seg[i + size] = vec[i + 1];
        for (int i = size - 1; i >= 1; --i)
            seg[i] = merge(seg[2 * i], seg[2 * i + 1]);
    }

    void update(int pos, int val) {
        pos += size - 1;
        seg[pos] = merge(seg[pos], val);
        for (pos >>= 1; pos; pos >>= 1)
            seg[pos] = merge(seg[2 * pos], seg[2 * pos + 1]);
    }

    int query(int l, int r) {
        l += size - 1;
        r += size - 1;
        int res = 0;
        while (l <= r) {
            if (l & 1) res = merge(res, seg[l++]);
            if (!(r & 1)) res = merge(res, seg[r--]);
            l >>= 1;
            r >>= 1;
        }
        return res;
    }
};

// Recursive Segment Tree
struct SegmentTree {
    using T = int;

    int n, neutral;
    vector<T> seg;

    T merge(T a, T b) { return min(a, b); }

    SegmentTree() {
    }

    SegmentTree(int _n) {
        n = _n, neutral = 0;
        seg.assign(4 * n + 5, neutral);
    }

    SegmentTree(const vector<T> &vec) {
        n = (int) vec.size() - 1, neutral = 0;
        seg.assign(4 * n + 5, neutral);
        build(1, 1, n, vec);
    }

    void build(int i, int l, int r, const vector<T> &vec) {
        if (l == r) {
            seg[i] = vec[l];
            return;
        }

        build(2 * i, l, (l + r) / 2, vec);
        build(2 * i + 1, (l + r) / 2 + 1, r, vec);
        seg[i] = merge(seg[2 * i], seg[2 * i + 1]);
    }

    void update(int i, int l, int r, int ind, int val) {
        if (ind < l || ind > r)
            return;

        if (l == r) {
            seg[i] = val;
            return;
        }

        update(2 * i, l, (l + r) / 2, ind, val);
        update(2 * i + 1, (l + r) / 2 + 1, r, ind, val);

        seg[i] = merge(seg[2 * i], seg[2 * i + 1]);
    }

    T query(int i, int l, int r, int lq, int rq) {
        if (r < lq || l > rq)
            return neutral;

        if (lq <= l && r <= rq)
            return seg[i];

        T a = query(2 * i, l, (l + r) / 2, lq, rq);
        T b = query(2 * i + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge(a, b);
    }
};

// ==========================================================================================

// Lazy Propagation, Lazy Segment Tree
struct LazySegmentTree {
    int n, neutral;
    vector<int> seg, lazy;

    LazySegmentTree() {}
    LazySegmentTree(int _n) {
        n = _n, neutral = 0;
        seg = lazy = vector<int>(4 * n + 5, neutral);
    }
    LazySegmentTree(const vector<int> &vec) {
        n = vec.size() - 1, neutral = 0;
        seg = lazy = vector<int>(4 * size + 2, neutral);
        build(1, 1, n, vec);
    }

    int merge(int a, int b) { return max(a, b); }
    void process(int i, int l, int r) {
        if (lazy[i] == 0)
            return;

        seg[n] += lazy[i];
        if (l != r)
            lazy[2 * i] += lazy[i], lazy[2 * i + 1] += lazy[i];
        lazy[i] = 0;
    }
    void build(int i, int l, int r, const vector<int> &vec) {
        process(i, l, r);

        if (l == r) {
            seg[i] = vec[l];
            return;
        }

        build(2 * i, l, (l + r) / 2, vec);
        build(2 * i + 1, (l + r) / 2 + 1, r, vec);
        seg[i] = merge(seg[2 * i], seg[2 * i + 1]);
    }
    void update(int i, int l, int r, int lq, int rq, int x) {
        process(i, l, r);

        if (r < lq || rq < l)
            return;

        if (lq <= l && r <= rq) {
            lazy[i] += x;
            process(i, l, r);
            return;
        }

        update(2 * i, l, (l + r) / 2, lq, rq, x);
        update(2 * i + 1, (l + r) / 2 + 1, r, lq, rq, x);

        seg[i] = merge(seg[2 * i], seg[2 * i + 1]);
    }
    int query(int i, int l, int r, int lq, int rq) {
        process(i, l, r);

        if (r < lq || l > rq)
            return neutral;

        if (lq <= l && r <= rq)
            return seg[i];

        int a = query(2 * i, l, (l + r) / 2, lq, rq);
        int b = query(2 * i + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge(a, b);
    }
};

// ==========================================================================================

// Maximum Sum Range Query
struct node {
    int prefix, middle, suffix, sum;

    node() { prefix = middle = suffix = sum = 0; }
    node(int a, int b, int c, int d) { prefix = a, middle = b, suffix = c, sum = d; }
};
struct MaximumSum {
    int n, INF = 1e18;
    vector<node> seg;

    MaximumSum() {}
    MaximumSum(int _n) {
        n = _n;
        seg = vector<node> (4 * n + 5, node());
    }

    node merge(node a, node b) {
        if (a.prefix == -INF)
            return b;
        if (b.prefix == -INF)
            return a;

        node ans(0, 0, 0, 0);
        ans.prefix = max({a.prefix, a.sum + b.prefix, a.sum + b.sum});
        ans.middle = max({a.middle, b.middle, a.suffix + b.prefix, a.sum + b.sum});
        ans.suffix = max({b.suffix, b.sum + a.suffix, a.sum + b.sum});
        ans.sum = a.sum + b.sum;

        return ans;
    }

    void build(int i, int l, int r, const vector<int>& vec) {
        if (l == r) {
            seg[i] = node(vec[l], vec[l], vec[l], vec[l]);
            return;
        }

        build(2 * i, l, (l + r) / 2, vec);
        build(2 * i + 1, (l + r) / 2 + 1, r, vec);
        seg[i] = merge(seg[2 * i], seg[2 * i + 1]);
    }
    void update(int i, int l, int r, int ind, int val) {
        if (ind < l || ind > r)
            return;

        if (l == r) {
            seg[i] = node(val, val, val, val);
            return;
        }

        update(2 * i, l, (l + r) / 2, ind, val);
        update(2 * i + 1, (l + r) / 2 + 1, r, ind, val);

        seg[i] = merge(seg[2 * i], seg[2 * i + 1]);
    }
    node query(int i, int l, int r, int lq, int rq) {
        if (l > rq || lq > r)
            return node(-INF, -INF, -INF, -INF);

        if (lq <= l && r <= rq)
            return seg[i];

        node a = query(2 * i, l, (l + r) / 2, lq, rq);
        node b = query(2 * i + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge(a, b);
    }
};

// ==========================================================================================

// Maximum Alternating Subsequence Sum
struct Node {
    int max_pp; // max of start with + and end with +
    int max_pn; // max of start with + and end with -
    int min_pp; // min of start with + and end with +
    int min_pn; // min of start with + and end with -

    Node() {}
    Node(int a, int b) {
        max_pp = a, min_pp = a;
        max_pn = b, min_pn = b;
    }
};
struct MaximumAlternatingSubsequenceSum {
    int n;
    vector<Node> seg;

    MaximumAlternatingSubsequenceSum(const vector<int> &vec) {
        int n = vec.size() - 1;
        seg = vector<Node>(4 * n + 5, Node(0, 0));
        build(1, 1, n, vec);
    }

    Node merge(Node a, Node b) {
        Node ans(0, 0);

        ans.max_pp = max({a.max_pp, b.max_pp, a.max_pp - b.min_pn, a.max_pn + b.max_pp});
        ans.max_pn = max({a.max_pn, b.max_pn, a.max_pn + b.max_pn, a.max_pp - b.min_pp});

        ans.min_pp = min({a.min_pp, b.min_pp, a.min_pp - b.max_pn, a.min_pn + b.min_pp});
        ans.min_pn = min({a.min_pn, b.min_pn, a.min_pn + b.min_pn, a.min_pp - b.max_pp});

        return ans;
    }

    void build(int i, int l, int r, const vector<int>& vec) {
        if (l == r) {
            seg[i] = Node(vec[l], 0);
            return;
        }

        build(2 * i, l, (l + r) / 2, vec);
        build(2 * i + 1, (l + r) / 2 + 1, r, vec);
        seg[i] = merge(seg[2 * i], seg[2 * i + 1]);
    }
    void update(int i, int l, int r, int ind, int val) {
        if (ind < l || ind > r)
            return;

        if (l == r) {
            seg[i] = Node(val, 0);
            return;
        }

        update(2 * i, l, (l + r) / 2, ind, val);
        update(2 * i + 1, (l + r) / 2 + 1, r, ind, val);

        seg[i] = merge(seg[2 * i], seg[2 * i + 1]);
    }
    Node query(int i, int l, int r, int lq, int rq) {
        if (l > rq || lq > r)
            return Node(0, 0);

        if (lq <= l && r <= rq)
            return seg[i];

        Node a = query(2 * i, l, (l + r) / 2, lq, rq);
        Node b = query(2 * i + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge(a, b);
    }
};

// ==========================================================================================

// Merge Sort Tree
struct MergeSortTree {
    int n;
    vector<vector<int> > seg;

    MergeSortTree() {
    }
    MergeSortTree(vector<int> &vec) {
        n = vec.size() - 1;
        seg = vector<vector<int>>(4 * n + 2, vector<int>{});
        build(1, 1, n, vec);
    }

    vector<int> merge(vector<int> &a, vector<int> &b) {
        vector<int> ans;
        int x = a.size(), y = b.size(), i = 0, j = 0;
        while (i < x && j < y) {
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

    void build(int i, int l, int r, const vector<int>& vec) {
        if (l == r) {
            seg[i] = vector<int>{vec[l]};
            return;
        }

        build(2 * i, l, (l + r) / 2, vec);
        build(2 * i + 1, (l + r) / 2 + 1, r, vec);
        seg[i] = merge(seg[2 * i], seg[2 * i + 1]);
    }
    int greater(int i, int l, int r, int lq, int rq, int val) {
        if (r < lq || l > rq)
            return 0;

        if (lq <= l && r <= rq) {
            int ind = upper_bound(seg[i].begin(), seg[i].end(), val) - seg[i].begin();
            return seg[i].size() - ind;
        }

        int a = greater(2 * i, l, (l + r) / 2, lq, rq, val);
        int b = greater(2 * i + 1, (l + r) / 2 + 1, r, lq, rq, val);

        return (a + b);
    }
    int less(int i, int l, int r, int lq, int rq, int val) {
        if (r < lq || l > rq)
            return 0;
        if (lq <= l && r <= rq) {
            int ind = lower_bound(seg[i].begin(), seg[i].end(), val) - seg[i].begin();
            return ind;
        }
        int a = less(2 * i, l, (l + r) / 2, lq, rq, val);
        int b = less(2 * i + 1, (l + r) / 2 + 1, r, lq, rq, val);
        return (a + b);
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
