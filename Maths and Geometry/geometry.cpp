#include <bits/stdc++.h>
#define int long long
#define double long double

using namespace std;

const double eps = 1e-9;
const double pi = acos(-1);
const double rotation = 360.0;
const double INF = 1e18;

#define X real()
#define Y imag()
#define point complex<double>

#define vec(a,b) ((b) - (a))
#define length(a) (hypot((a).X, (a).Y))
#define normalize(a) ((a)/length(a))

#define dot(a,b) ((conj(a)*(b)).X)
#define cross(a,b) ((conj(a)*(b)).Y)
#define same(a,b) (dot(vec(a,b),vec(b,a)) < eps)
#define square(a) (dot((a),(a)))

#define rotate0(a,ang) ((a)*exp(point(0,ang)))
#define rotate_about(a,ang,about) (rotate0(vec(about,a),ang)+about)

#define reflect0(a,b) (conj((a)/(b))*(b))
#define reflect_about(p, p0, p1) (reflect0(vec(p0, p), vec(p0, p1)) + p0)

double to_radians(double degree) { return degree / 180.0 * pi; }
double to_degree(double radian) { return radian / pi * 180.0; }
int dcmp(double x, double y) { return (fabs(x - y) < eps) ? 0 : (x < y ? -1 : 1); }

double angle(point a) { return (atan2((a).Y, (a).X)); }
double angle(point v, point w) { return acos(min((double)1, max((double)-1, dot(v, w) / abs(v) / abs(w)))); } // angle between v, w
point translate(point v, point p) { return v + p; } // translates v by p
point scale(point v, double factor, point center = {0, 0}) { return center + (v - center) * factor; } // scales v by factor around center
point linear_transformation(point p, point fp, point q, point fq, point r) { return fp + (r - p) * ((fq - fp) / (q - p)); } // p->fp, q->fq, r->?
point perp(point p) { return point(-p.Y, p.X); }

int orient(point a, point b, point c) { return dcmp(cross(b - a, c - a), 0); } // a->b->c to left or right
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

bool is_collinear(point a, point b, point c) {
    return dcmp(fabs(cross(b - a, c - a)), 0) == 0;
}
bool is_point_on_ray(point a, point b, point c) // is point c on ray ab
{
    if (!is_collinear(a, b, c))
        return false;
    return dcmp(dot(b - a, c - a), 0) >= 0;

    // if (length(a - c) == 0)
    //     return true;
    // return same(normalize(b - a), normalize(c - a));
}
bool is_point_on_segment(point a, point b, point c) // is point c on segment ab
{
    return is_point_on_ray(a, b, c) && is_point_on_ray(b, a, c);

    // double ab = length(b - a), ac = length(c - a), bc = length(c - b);
    // return dcmp(ab, ac + bc) == 0;
}
double distance_to_line(point a, point b, point c) // distance between point c and line ab
{
    double d = cross(b - a, c - a) / length(a - b);
    return fabs(d);
}
double distance_to_segment(point a, point b, point c) // distance between point c and segment ab
{
    if (dcmp(dot(b - a, c - a), 0) == -1)
        return length(c - a);
    if (dcmp(dot(a - b, c - b), 0) == -1)
        return length(c - b);
    return distance_to_line(a, b, c);
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

pair<double, point > circle(point a, point b, point c) {
    if (is_collinear(a, b, c))
        return {-1, point(0, 0)};

    point m1 = (a + b) / (double) 2.0, m2 = (a + c) / (double) 2.0;
    point v1 = b - a, v2 = c - a;
    point p1(-v1.Y, v1.X), p2(-v2.Y, v2.X);
    point end1 = m1 + p1, end2 = m2 + p2;
    point center;
    bool ok = intersect(line(m1, end1), line(m2, end2), center);
    return {length(a - center), center};
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

double triangle_area(point a, point b, point c) { return abs(cross(b - a, c - a)) / 2; } // calc area of triangle abc
double polygon_area(vector<point>& p) { // calc area of a polygon
    int n = p.size();
    double area = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += cross(p[i], p[j]);
    }
    return abs(area) / 2;
}
bool above(point a, point p) { return dcmp(p.Y, a.Y) >= 0; } // is point p in the above half plane of a
bool ray_cross(point a, point p, point q) { // does segment pq crosses the + horizontal ray of a
    return ((int)above(a, q) - (int)above(a, p)) * orient(a, p, q) > 0;
}
bool point_in_polygon(vector<point>& p, point a, bool strict = true) { // is point inside a polygon O(n)
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
bool is_polygon_simple(vector<point > &p) {
    bool ok = true;
    p.push_back(p[0]);
    for (int i = 0; i < (int) p.size() - 1; i++) {
        for (int j = i + 2; j < (int) p.size() - 1; j++) {
            auto ret = improper_segments_intersection(p[i], p[i + 1], p[j], p[j + 1]);
            if (ret.size() && !(i == 0 && j == (int) p.size() - 2))
                ok = false;
        }
    }
    p.pop_back();
    return ok;
}
bool is_polygon_convex(const vector<point>& points) { // is polygon convex (co-linear lines are considered convex)
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
point centroid(vector<point > &p) {
    p.push_back(p[0]);

    double area = 0, x = 0, y = 0;
    for (int i = 0; i < (int) p.size() - 1; i++) {
        double c = cross(p[i], p[i + 1]);
        area += c;
        x += (p[i].X + p[i + 1].X) * c;
        y += (p[i].Y + p[i + 1].Y) * c;
    }
    p.pop_back();
    area /= 2, x /= 6 * area, y /= 6 * area;

    if (area == 0)
        return (p[0] + p.back()) * (double) 0.5;

    if (dcmp(x, 0) == 0)
        x = 0;
    if (dcmp(y, 0) == 0)
        y = 0;

    return point(x, y);
}
pair<vector<point >, vector<point > > polygon_cut(vector<point > &p, point a, point b) {
    p.push_back(p[0]);
    vector<point > left, right;

    for (int i = 0; i < (int) p.size(); i++) {
        if (dcmp(cross(b - a, p[i] - a), 0) >= 0)
            right.push_back(p[i]);

        auto ret = improper_segments_intersection(a, b, p[i], p[i + 1]);
        if (ret.size()) {
            point j;
            bool ok = intersect(line(a, b), line(p[i], p[i + 1]), j);
            right.push_back(j);
            left.push_back(j);
        }

        if (dcmp(cross(b - a, p[i] - a), 0) <= 0)
            left.push_back(p[i]);
    }
    return {left, right};
}

// Is point in a polygon (clockwise convex polygon)
int is_point_inside_polygon(const point& p, const vector<point>& poly) {
    int min_x = 1e10, max_x = -1e10;
    int min_y = 1e10, max_y = -1e10;

    for (auto& po : poly) {
        int x = round(po.X);
        int y = round(po.Y);

        min_x = min(min_x, x);
        max_x = max(max_x, x);
        min_y = min(min_y, y);
        max_y = max(max_y, y);
    }

    int x = round(p.X);
    int y = round(p.Y);

    if (x < min_x || x > max_x || y < min_y || y > max_y)
        return 0;

    int n = poly.size();
    int l = 1, r = n - 2, ans = -1;

    while (l <= r) {
        int mid = (l + r) / 2;
        int j = (mid + 1) % n;
        if (orient(poly[0], poly[mid], p) <= 0 && orient(poly[0], poly[j], p) >= 0) {
            ans = mid;
            break;
        }
        if (orient(poly[0], poly[mid], p) >= 0 && orient(poly[0], poly[j], p) <= 0) {
            ans = mid;
            break;
        }
        if (orient(poly[0], poly[mid], p) >= 0 && orient(poly[0], poly[j], p) >= 0) {
            r = mid - 1;
        }
        else if (orient(poly[0], poly[mid], p) <= 0 && orient(poly[0], poly[j], p) <= 0) {
            l = mid + 1;
        }
    }

    if (ans == -1)
        return 0;

    int j = (ans + 1) % n;
    if (orient(poly[ans], poly[j], p) > 0)
        return 0;
    return 1;
}

// Closest Pair
double closest_pair(vector<point>& p, pair<int, int>& ret) {
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
    bool operator<(const segment &b) const {
        double x = max(min(p.X, q.X), min(p.X, b.q.X));
        return dcmp(get_y(x), b.get_y(x)) == -1;
    }
};
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

set<segment> s;
vector<set<segment>::iterator> where;
set<segment>::iterator prev(set<segment>::iterator it) { return it == s.begin() ? s.end() : --it; }
set<segment>::iterator next(set<segment>::iterator it) { return ++it; }

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
pair<int, int> intersecting_segments(const vector<segment> &segments) {
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

// Convex Hull (clockwise)
bool cw(point a, point b, point c, bool include_collinear) {
    int o = orient(a, b, c);
    return o < 0 || (include_collinear && o == 0);
}
bool collinear(point a, point b, point c) { return orient(a, b, c) == 0; }
void convex_hull(vector<point>& a, bool include_collinear = false) {
    point p0 = *min_element(a.begin(), a.end(), [](point a, point b) {
        return make_pair(a.Y, a.X) < make_pair(b.Y, b.X);
    });
    sort(a.begin(), a.end(), [&p0](const point& a, const point& b) {
        int o = orient(p0, a, b);
        if (o == 0)
            return  dcmp(dot(p0 - a, p0 - a), dot(p0 - b, p0 - b)) < 0;
        return o < 0;
    });
    if (include_collinear) {
        int i = (int)a.size()-1;
        while (i >= 0 && collinear(p0, a[i], a.back())) i--;
        reverse(a.begin()+i+1, a.end());
    }

    vector<point> st;
    for (int i = 0; i < (int)a.size(); i++) {
        while (st.size() > 1 && !cw(st[st.size()-2], st.back(), a[i], include_collinear))
            st.pop_back();
        st.push_back(a[i]);
    }

    if (include_collinear == false && st.size() == 2 && st[0] == st[1])
        st.pop_back();

    a = st;
}

// Convex Hull (counterclockwise, another implementation)
bool comp(const point& a, const point& b) {
    if (dcmp(a.X, b.X) != 0)
        return dcmp(a.X, b.X) < 0;
    return dcmp(a.Y, b.Y) < 0;
}
vector<point> convex_hull(vector<point>& p) {
    sort(p.begin(), p.end(), comp);

    int n = p.size();
    point p1 = p[0], p2 = p.back();

    vector<point> up = {p1}, down = {p1};
    for (int i = 1; i < n; i++) {
        if (i == n - 1 || orient(p1, p[i], p2) <= 0) {
            while (up.size() > 1 && orient(up[up.size() - 2], up.back(), p[i]) >= 0)
                up.pop_back();
            up.push_back(p[i]);
        }
        if (i == n - 1 || orient(p1, p[i], p2) >= 0) {
            while (down.size() > 1 && orient(down[down.size() - 2], down.back(), p[i]) <= 0)
                down.pop_back();
            down.push_back(p[i]);
        }
    }

    return up;
    vector<point> res;
    for (int i = 0; i < down.size(); i++)
        res.push_back(down[i]);
    for (int i = (int)up.size() - 2; i >= 1; i--)
        res.push_back(up[i]);
    return res;
}

// Rectangles Union
int rectangles_union(const vector<array<int, 4>>& recs) {
    const int N = 3e4 + 5;
    struct event {
        int x, y1, y2, type;
        event() {}
        event(int _x, int _y1, int _y2, int _type) : x(_x), y1(_y1), y2(_y2), type(_type) {}

        bool operator<(const event& e) const {
            if (x != e.x)
                return x < e.x;
            return type > e.type;
        }
    };
    struct LazySegmentTree {
        int n, neutral;
        vector<int> seg, full;

        LazySegmentTree() {}
        LazySegmentTree(int _n) {
            n = _n, neutral = 0;
            seg = full = vector<int>(4 * n + 5, neutral);
        }

        int get(int i, int l, int r) {
            if (full[i])
                return r - l + 1;
            return seg[i];
        }
        void update(int i, int l, int r, int lq, int rq, int x) {
            if (r < lq || rq < l)
                return;

            if (lq <= l && r <= rq) {
                full[i] += x;
                return;
            }

            int mid = (l + r) / 2;
            update(2 * i, l, mid, lq, rq, x);
            update(2 * i + 1, mid + 1, r, lq, rq, x);

            seg[i] = get(2 * i, l, mid) + get(2 * i + 1, mid + 1, r);
        }
        int query(int i, int l, int r, int lq, int rq) {
            if (r < lq || l > rq)
                return 0;

            if (lq <= l && r <= rq)
                return seg[i];

            int mid = (l + r) / 2;
            int a = query(2 * i, l, mid, lq, rq);
            int b = query(2 * i + 1, mid + 1, r, lq, rq);

            return a + b;
        }
    };

    int n = recs.size();
    vector<event> v;
    for (int i = 0; i < n; i++) {
        auto [x1, y1, x2, y2] = recs[i];
        v.push_back(event(x1, y1, y2 - 1, +1));
        v.push_back(event(x2, y1, y2 - 1, -1));
    }
    sort(v.begin(), v.end());

    int ans = 0;
    int last = v[0].x;
    LazySegmentTree tree(N);

    for (auto& c : v) {
        ans += (c.x - last) * tree.query(1, 1, N, 1, N);
        tree.update(1, 1, N, c.y1, c.y2, c.type);
        last = c.x;
    }

    return ans;
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);



    return 0;
}