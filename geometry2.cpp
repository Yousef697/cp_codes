#include <bits/stdc++.h>
#define int long long
#define double long double

using namespace std;

typedef complex<double> point;
#define X real()
#define Y imag()

const double eps = 1e-9;
const double pi = acos(-1);
const double INF = 1e18;

int dcmp(double x, double y) { return (fabs(x - y) < eps ? 0 : x > y ? 1 : -1); }

point translate(point v, point p) { return v + p; }// translates v by p
point scale(point v, double factor, point center = {0, 0}) { return center + (v - center) * factor; } // scales v by factor around center
point rotate0(point v, double angle) { return point( v.X * cos(angle) - v.Y * sin(angle), v.X * sin(angle) + v.Y * cos(angle)); } // rotate v by angle around origin
point rotate(point v, double angle, point about) { return about + rotate0(v - about, angle); } // rotate v by angle around about
point linear_transformation(point p, point fp, point q, point fq, point r) { return fp + (r - p) * ((fq - fp) / (q - p)); } // p->fp, q->fq, r->?
point perp(point p) { return point(-p.Y, p.X); }

double dot(point v, point w) { return (conj(v) * w).X; } // dot product of v, w
double cross(point v, point w) { return (conj(v) * w).Y; } // cross product of v, w
int orient(point a, point b, point c) { return dcmp(cross(b - a, c - a), 0); } // a->b->c to left or right
double angle(point v, point w) { return acos(min((double)1, max((double)-1, dot(v, w) / abs(v) / abs(w)))); } // angle between v, w
bool in_angle(point a, point b, point c, point p) { // is point p inside the angle bac
    assert(orient(a, b, c) != 0);
    if (orient(a, b, c) == -1)
        swap(b, c);
    return orient(a, b, p) >= 0 && orient(a, c, p) <= 0;
}
double oriented_angle(point a, point b, point c) { // angle bac
    if (orient(a, b, c) >= 0)
        return angle(b - a, c - a);
    return 2 * pi - angle(b - a, c - a);
}
bool is_convex(const vector<point>& points) { // is polygon convex (co-linear lines are considered convex)
    int n = points.size();
    bool positive = false, negative = false;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        int k = (j + 1) % n;

        int ret = orient(points[i], points[j], points[k]);
        if (ret < 0)
            negative = true;
        else if (ret > 0)
            positive = true;
    }
    return !(positive && negative);
}

point around = {1, 1}, about = {1, 1};
bool half(point p) {
    assert(p.X != 0 || p.Y != 0);
    return dcmp(p.Y, 0) == 1 || (dcmp(p.Y, 0) == 0 && dcmp(p.X, 0) == -1);
}
bool half2(point p) {
    return dcmp(cross(about, p), 0) == -1 ||
        (dcmp(cross(about, p), 0) == 0 && dcmp(dot(about, p), 0) == -1);
}
void polar_sort(vector<point>& points) {
    // polar sort around origin
    sort(points.begin(), points.end(), [](point& v, point& w) {
        return make_tuple(half(v), (double)0) < make_tuple(half(w), cross(v, w));
    });

    // polar & length sort around origin
    sort(points.begin(), points.end(), [](point& v, point& w) {
        return make_tuple(half(v), (double)0, dot(v, v)) < make_tuple(half(w), cross(v, w), dot(w, w));
    });

    // polar sort around a point
    sort(points.begin(), points.end(), [](point& v, point& w) {
        return make_tuple(half(v - around), (double)0) < make_tuple(half(w - around), cross(v - around, w - around));
    });

    // polar sort with respect to about
    sort(points.begin(), points.end(), [](point& v, point& w) {
        return make_tuple(half2(v), (double)0) < make_tuple(half2(w), cross(v, w));
    });
}

struct line { // ax + by = c
    point v;
    double c;

    line(point v, double c) : v(v), c(c) {}
    line(double a, double b, double c) : v({b, -a}), c(c) {}
    line(point p, point q) : v(q - p), c(cross(v, p)) {}

    double side(point p) { return cross(v, p) - c; } // point on which side of the line (+ = positive, - = negative)
    double dist(point p) { return abs(side(p)) / abs(v); } // distance between point and line
    line perp_through(point p) { return {p, p + perp(v)}; } // perpendicular line through the point p
    bool compare_proj(point p, point q) { return dot(v, p) < dot(v, q); } // p comes before q if it appears on the line before q
    line translate(point t) { return {v, c + cross(v, t)}; } // translate the line in the direction of t
    line shift_left(double dist) { return {v, c + dist * abs(v)}; } // translates the line with distance dist to left
    point projection(point p) { return p - perp(v) * side(p) / dot(v, v); } // projection of a point on the line
    point reflect(point p) { return p - (double)2 * perp(v) * side(p) / dot(v, v); } // reflection of a point with respect to the line
};
bool intersect(line l1, line l2, point& ret) { // intersection of two lines
    double d = cross(l1.v, l2.v);
    if (dcmp(d, 0) == 0)
        return false;
    ret = (l1.c * l2.v - l2.c * l1.v) / d;
    return true;
}
line bisector(line l1, line l2, bool interior) { // bisector of two lines
    assert(dcmp(cross(l1.v, l2.v), 0) != 0);
    double sign = interior ? 1 : -1;
    return line(l2.v / abs(l2.v) + l1.v / abs(l1.v) * sign, l2.c / abs(l2.v) + l1.c / abs(l1.v) * sign);
}

bool in_circle(point a, point b, point p) { return dcmp(dot(a - p, b - p), 0) <= 0; } // is point p inside the circle with diameter ab
bool point_on_segment(point a, point b, point p) { return orient(a, b, p) == 0 && in_circle(a, b, p); } // is point p on segment ab
bool proper_segments_intersection(point a, point b, point c, point d, point& ret) { // are two segments intersect in a non-endpoint point
    int oa = orient(c, d, a);
    int ob = orient(c, d, b);
    int oc = orient(a, b, c);
    int od = orient(a, b, d);

    if (oa * ob < 0 && oc * od < 0) {
        double x = cross(d - c, a - c);
        double y = cross(d - c, b - c);
        ret = (a * x - b * y) / (y - x);
        return true;
    }
    return false;
}
struct cmp {
    bool operator()(const point& a, const point& b) const {
        return dcmp(a.X, b.X) == -1 || dcmp(a.X, b.X) == 0 && dcmp(a.Y, b.Y) == -1;
    }
};
set<point, cmp> improper_segments_intersection(point a, point b, point c, point d) { // intersection of two segments if they do not intersect properly
    point ret;
    if (proper_segments_intersection(a, b, c, d, ret)) { return {ret}; }

    set<point, cmp> intersections;
    if (point_on_segment(c, d, a)) intersections.insert(a);
    if (point_on_segment(c, d, b)) intersections.insert(b);
    if (point_on_segment(a, b, c)) intersections.insert(c);
    if (point_on_segment(a, b, d)) intersections.insert(d);
    return intersections;
};
double point_segment_distance(point a, point b, point p) { // minimum distance between segment and a point
    if (!(dcmp(a.X, b.X) == 0 && dcmp(a.Y, b.Y) == 0)) {
        line l(a, b);
        if (l.compare_proj(a, p) && l.compare_proj(b, p)) {
            return l.dist(p);
        }
    }
    return min(abs(p - a), abs(p - b));
}
double segment_segment_distance(point a, point b, point c, point d) { // minimum distance between segment and a segment
    point x;
    if (proper_segments_intersection(a, b, c, d, x)) {
        return 0;
    }
    double mn = point_segment_distance(a, b, c);
    mn = min(mn, point_segment_distance(a, b, d));
    mn = min(mn, point_segment_distance(c, d, a));
    mn = min(mn, point_segment_distance(c, d, b));
    return mn;
}

double triangle_area(point a, point b, point c) { return cross(b - a, c - a) / 2; } // calc area of triangle abc
double polygon_area(vector<point>& p) { // calc area of a polygon
    int n = p.size();
    double area = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += cross(p[i], p[j]);
    }
    return abs(area) / 2;
}
bool above(point a, point p) { return dcmp(p.Y, a.Y) >= 0; } // is point p in the above half of a
bool ray_cross(point a, point p, point q) { // does segment pq crosses the + horizontal ray of a
    return ((int)above(a, q) - (int)above(a, p)) * orient(a, p, q) > 0;
}
bool point_in_polygon(vector<point>& p, point a, bool strict = true) { // is point inside a polygon
    int n = p.size(), crosses = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        if (point_on_segment(p[i], p[j], a))
            return !strict;
        crosses += ray_cross(a, p[i], p[j]);
    }
    return crosses & 1;
}
double angle_travelled(point a, point p, point q) { // calc angle paq
    double ang = angle(p - a, q - a);
    if (orient(a, p, q))
        return ang;
    return -ang;
}
int winding_number(vector<point>& p, point a) { // calc winding number that test if a point inside a polygon
    int n = p.size();
    double ang = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        ang += angle_travelled(a, p[i], p[j]);
    }
    return round(ang / (2 * pi));
}

point circumcircle(point a, point b, point c) { // calc the center of the circle that passes through a, b, c
    point ab = b - a;
    point ac = c - a;

    assert(dcmp(cross(ab, ac), 0) != 0);
    return a + perp(ab * dot(ac, ac) - ac * dot(ab, ab)) / (2.0 * cross(ab, ac));
}
int circle_line_intersection(point center, double r, line l, pair<point, point>& ret) { // calc point of intersection of a circle and a line
    double h2 = r * r - l.dist(center) * l.dist(center);
    if (dcmp(h2, 0) >= 0) {
        point p = l.projection(center);
        point h = l.v * sqrt(h2) / abs(l.v);
        ret = {p - h, p + h};
    }
    return 1 + dcmp(h2, 0);
}
int circle_circle_intersection(point cen1, double r1, point cen2, double r2, pair<point, point>& ret) {  // calc point of intersection of a circle and a circle
    point d = cen2 - cen1;
    double d2 = dot(d, d);

    if (dcmp(d2, 0) == 0) {
        assert(dcmp(r1, r2) != 0);
        return 0;
    }

    double pd = (d2 + r1 * r1 - r2 * r2) / 2;
    double h2 = r1 * r1 - pd * pd / d2;

    if (dcmp(h2, 0) >= 0) {
        point p = cen1 + d * pd / d2;
        point h = perp(d) * sqrt(h2 / d2);
        ret = {p - h, p + h};
    }
    return 1 + dcmp(h2, 0);
}
int circles_tangents(point cen1, double r1, point cen2, double r2, bool inner, vector<pair<point, point>>& ret) { // calc the common tangents of two circles
    if (inner)
        r2 = -r2;
    point d = cen2 - cen1;
    double dr = r2 - r1, d2 = dot(d, d), h2 = d2 - dr * dr;
    if (dcmp(d2, 0) == 0 || dcmp(h2, 0) == -1) {
        assert(dcmp(h2, 0) != 0);
        return 0;
    }
    for (double sign : {-1, 1}) {
        point v = (d * dr + perp(d) * sqrt(h2) * sign) / d2;
        ret.push_back({cen1 + v * r1, cen2 + v * r2});
    }
    return 1 + dcmp(h2, 0);
}

point centroid(vector<point>& p) {

    int n = p.size();
    double area = 0, x = 0, y = 0;

    p.push_back(p[0]);
    for (int i = 0; i < n; i++) {
        double c = cross(p[i], p[i + 1]);
        area += c;
        x += (p[i].X + p[i + 1].X) * c;
        y += (p[i].Y + p[i + 1].Y) * c;
    }
    p.pop_back();
    area /= 2, x /= 6 * area, y /= 6 * area;

    if (area == 0)
        return (p[0] + p.back()) * (double)0.5;

    if (dcmp(x, 0) == 0)
        x = 0;
    if (dcmp(y, 0) == 0)
        y = 0;

    return point(x, y);
}

// Closest Pair
struct compare_x {
    bool operator() (const point& a, const  point& b) const {
        if (dcmp(a.X, b.X) != 0)
            return dcmp(a.X, b.X) == -1;
        return dcmp(a.Y, b.Y) == -1;
    }
};
struct compare_y {
    bool operator() (const point& a, const  point& b) const {
        if (dcmp(a.Y, b.Y) != 0)
            return dcmp(a.Y, b.Y) == -1;
        return dcmp(a.X, b.X) == -1;
    }
};
double closest_pair(vector<point>& p, pair<int, int>& ret) {
    int n = p.size();
    double ans = INF;
    vector<pair<int, int>> q;
    for (int i = 0; i < n; i++) {
        int x = round(p[i].X);
        int y = round(p[i].Y);
        q.push_back({x, y});
    }
    map<pair<int, int>, int> freq;
    for (int i = 0; i < n; i++)
        freq[q[i]] = i;

    int left = 0;
    point p1(-1, -1), p2(-1, -1);
    multiset<point, compare_y> active;
    sort(p.begin(), p.end(), compare_x());

    for (int right = 0; right < n; right++) {
        while (!active.empty() && left < right && dcmp(p[right].X - p[left].X, ans) == 1) {
            active.erase(active.find(p[left])), left++;
        }

        auto st = active.lower_bound(point(-INF, p[right].Y - ans));
        auto en = active.upper_bound(point(-INF, p[right].Y + ans));

        for (; st != en; st++) {
            point vec = p[right] - *st;
            double cand = dot(vec, vec);
            cand = sqrt(cand);
            if (dcmp(cand, ans) == -1) {
                ans = cand;
                p1 = *st;
                p2 = p[right];
            }
        }
        active.insert(p[right]);
    }

    int x = round(p1.X), y = round(p1.Y);
    ret.first = freq[{x, y}];
    x = round(p2.X), y = round(p2.Y);
    ret.second = freq[{x, y}];

    return ans;
}

// Intersecting Segments
struct segment {
    point p, q;
    int id;

    segment() {}

    segment(point p, point q, int id) : p(p), q(q), id(id) {}

    double get_y(double x) const {
        if (dcmp(p.X, q.X) == 0)
            return p.Y;
        return p.Y + (q.Y - p.Y) * (x - p.X) / (q.X - p.X);
    }
};
bool operator<(const segment &a, const segment &b) {
    double x = max(min(a.p.X, a.q.X), min(b.p.X, b.q.X));
    return dcmp(a.get_y(x), b.get_y(x)) == -1;
}
struct event {
    double x;
    int tp, id;

    event() {}

    event(double x, int tp, int id) : x(x), tp(tp), id(id) {}

    bool operator<(const event &e) const {
        if (dcmp(x, e.x) != 0)
            return x < e.x;
        return tp > e.tp;
    }
};

bool intersect1d(double l1, double r1, double l2, double r2) {
    if (l1 > r1)
        swap(l1, r1);
    if (l2 > r2)
        swap(l2, r2);
    return dcmp(max(l1, l2), min(r1, r2)) <= 0;
}
bool intersect(const segment &a, const segment &b) {
    return intersect1d(a.p.X, a.q.X, b.p.X, b.q.X) &&
           intersect1d(a.p.Y, a.q.Y, b.p.Y, b.q.Y) &&
           orient(a.p, a.q, b.p) * orient(a.p, a.q, b.q) <= 0 &&
           orient(b.p, b.q, a.p) * orient(b.p, b.q, a.q) <= 0;
}

set<segment> s;
vector<set<segment>::iterator> where;
set<segment>::iterator prev(set<segment>::iterator it) { return it == s.begin() ? s.end() : --it; }
set<segment>::iterator next(set<segment>::iterator it) { return ++it; }

pair<int, int> solve(const vector<segment> &segments) {
    int n = (int) segments.size();
    vector<event> e;
    for (int i = 0; i < n; ++i) {
        e.push_back(event(min(segments[i].p.X, segments[i].q.X), +1, i));
        e.push_back(event(max(segments[i].p.X, segments[i].q.X), -1, i));
    }
    sort(e.begin(), e.end());

    s.clear();
    where.resize(segments.size());

    for (size_t i = 0; i < e.size(); ++i) {
        int id = e[i].id;
        if (e[i].tp == +1) {
            set<segment>::iterator nxt = s.lower_bound(segments[id]), prv = prev(nxt);
            if (nxt != s.end() && intersect(*nxt, segments[id]))
                return make_pair(nxt->id, id);
            if (prv != s.end() && intersect(*prv, segments[id]))
                return make_pair(prv->id, id);
            where[id] = s.insert(nxt, segments[id]);
        } else {
            set<segment>::iterator nxt = next(where[id]), prv = prev(where[id]);
            if (nxt != s.end() && prv != s.end() && intersect(*nxt, *prv))
                return make_pair(prv->id, nxt->id);
            s.erase(where[id]);
        }
    }

    return make_pair(-1, -1);
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    
    
    return 0;
}
