#include <bits/stdc++.h>
#include <complex>
#include <cmath>

using namespace std;

const double eps = 1e-9;
const double pi = acos(-1);
const double rotation = 360.0;

#define X real()
#define Y imag()
#define point complex<double>

#define angle(a) (atan2((a).Y, (a).X))
#define vec(a,b) ((b) - (a))
#define length(a) (hypot((a).X, (a).Y))
#define normalize(a) ((a)/length(a))

#define dot(a,b) ((conj(a)*(b)).X)
#define cross(a,b) ((conj(a)*(b)).Y)
#define same(a,b) (dot(vec(a,b),vec(b,a)) < eps)
#define square(a) (dot(a,a))

#define rotate(a,ang) ((p)*exp(point(0,ang)))
#define rotate(a,ang,about) (rotate(vec(about,a),ang)+about)

#define reflect(a,b) (conj(a/b)*b)

double to_radians(double degree) {
    return degree / 180.0 * pi;
}
double to_degree(double radian) {
    return radian / pi * 180.0;
}

int dcmp(double x, double y) {
    if (fabs(x - y) < eps) return 0;
    return (x < y ? -1 : 1);
}

bool is_collinear(point a, point b, point c) {
    return dcmp(cross(b - a, c - a), 0) == 0;
}
bool is_point_on_ray(point a, point b, point c) // is point c on ray ab
{
    if (!is_collinear(a, b, c))
        return false;
    return dcmp(dot(b - a, c - a), 0) == 1;

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
    double d = cross(a - c, b - c) / length(b - a);
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
point lines_intersect(point a, point b, point c, point d) {
    double d1 = cross(a - b, d - c);
    double d2 = cross(a - c, d - c);
    double d3 = cross(a - b, a - c);

    if (fabs(d1) < eps)
        return point(INT_MIN, INT_MIN);

    double t1 = d2 / d1;
    return a + (b - a) * t1;
}
int counter_clockwise(point a, point b, point c) {
    point v1 = b - a, v2 = c - b, v3 = c - a;

    if (cross(v1, v2) > eps)
        return 1;
    if (cross(v1, v2) < -eps)
        return -1;

    if (is_point_on_segment(a, b, c))
        return 0;

    if (same(normalize(v1), normalize(v2)))
        return 1;
    return -1;
}
bool are_segments_intersect(point a, point b, point c, point d) {
    // bool x = same(a, b), y = same(c, d);
    // if (x && y) return same(a, c);
    // if (x) return counter_clockwise(c, d, a);
    // if (y) return counter_clockwise(a, b, c);
    //
    // x = counter_clockwise(a, b, c) * counter_clockwise(a, b, d) == -1;
    // y = counter_clockwise(c, d, a) * counter_clockwise(c, d, b) == -1;
    // return x && y;
    point intersect = lines_intersect(a, b, c, d);
    return is_point_on_segment(a, b, intersect) && is_point_on_segment(c, d, intersect);
}

pair<double, point> circle(point a, point b, point c) {
    if (is_collinear(a, b, c))
        return {-1, point(0, 0)};

    point m1 = (a + b) / 2.0, m2 = (a + c) / 2.0;
    point v1 = b - a, v2 = c - a;
    point p1(-v1.Y, v1.X), p2(-v2.Y, v2.X);
    point end1 = m1 + p1, end2 = m2 + p2;
    point center = lines_intersect(m1, end1, m2, end2);
    return {length(a - center), center};
}
vector<point> circle_line_intersection(double r, point c, point a, point b) // radius, center, line ab
{
    double A = dot(b - a, b - a);
    double B = 2 * dot(b - a, a - c);
    double C = dot(a - c, a - c) - r * r;
    double disc = B * B - 4 * A * C;

    vector<point> res;
    if (dcmp(disc, 0) >= 0) {
        if (dcmp(disc, 0) == 0)
            disc = 0;
        double t1 = (-B + sqrt(disc)) / (2.0 * A);
        double t2 = (-B - sqrt(disc)) / (2.0 * A);

        res.push_back(a + t1 * (b - a));
        if (dcmp(t1, t2) != 0)
            res.push_back(a + t2 * (b - a));
    }
    return res;
}
vector<point> circle_circle_intersection(double r1, point c1, double r2, point c2) {
    // same circle with positive radius
    if (same(c1, c2) && dcmp(r1, r1) == 0 && dcmp(r1, 0) == 1)
        return vector<point>(3, c1);

    auto get_angle = [&](double a, double b, double c) {
        return acos((b*b + c*c - a*a) / (2 * b * c));
    };

    double ang1 = angle(c2 - c1), ang2 = get_angle(r2, r1, length(c2 - c1));
    if (isnan(ang2))
        ang2 = 0;

    vector<point> res;
    point p1 = c1 + polar(r1, ang1 + ang2);

    if (dcmp(dot(p1 - c1, p1 - c1), r1 * r1) != 0 || dcmp(dot(p1 - c2, p1 - c2), r2 * r2) != 0)
        return vector<point>();

    res.push_back(p1);
    point p2 = c1 + polar(r1, ang1 - ang2);
    if (!same(p1, p2))
        res.push_back(p2);

    return res;
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    complex<double> num1(2, 3);
    cout << num1 << "\n"; // (2,3)
    cout << num1.real() << " " << num1.imag() << "\n"; // 2 3

    complex<double> num2(1, 1);
    cout << "Abs: " << abs(num2) << "\n";
    cout << "Arg(rad): " << arg(num2) << "\n";
    cout << "Arg(rad): " << to_degree(arg(num2)) << "\n";
    cout << "Norm: " << norm(num2) << "\n----------\n";

    complex<double> num3 = polar(1.41421, 0.785398);
    cout << num3 << "\n----------\n";

    complex<double> zero;
    // complex<double> x_part = 10; does not work?!

    complex<double> a(1, 2), b(3, 4);
    cout << a << " + " << b << " = " << a + b << "\n";
    cout << a << " - " << b << " = " << a - b << "\n";
    cout << a << " * " << b << " = " << a * b << "\n";
    cout << a << " * " << 2.0 << " = " << a * 2.0 << "\n";
    cout << a << " / " << 2.0 << " = " << a / 2.0 << "\n----------\n";

    complex<double> i(0, 1);
    complex<double> i1(-1, 0);
    complex<double> i2(2, 3);
    cout << sqrt(i1) << "\n";
    cout << conj(i2) << "\n";
    cout << pow(i2, 2) << "\n";
    cout << exp(i * pi) << "\n----------\n";

    point A(2, 1);
    point B(6, 2);
    point C(4, 3); // (8, 3)
    point D(5, -2);
    cout << lines_intersect(A, B, C, D) << "\n----------\n";
    cout << are_segments_intersect(A, B, C, D) << "\n----------\n";

    point q(2, 0);
    point w(0, 2);
    point e(4, 4);
    auto p = circle(q, w, e);
    cout << p.first << ", " << p.second << "\n----------\n";

    auto ret1 = circle_line_intersection(3, point(2, 1), point(2, 3), point(-2, -2));
    auto ret2 = circle_line_intersection(3, point(0, 0), point(-4.24264068712, 0), point(0, 4.24264068712));
    auto ret3 = circle_line_intersection(3, point(0, 0), point(-4, 0), point(0, 8));
    for (auto i : ret1)
        cout << i << "\n";
    cout << "-\n";
    for (auto i : ret2)
        cout << i << "\n";
    cout << "-\n";
    for (auto i : ret3)
        cout << i << "\n";
    cout << "----------\n";

    return 0;
}
