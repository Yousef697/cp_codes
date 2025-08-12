/*
    Intro:
        we can solve triangles with this two rules:
        rule of sins: sin(A) / a = sin(B) / b = sin(C) / c
        rule of cosines: a^2 = b^2 + c^2 - 2bc * cos(A)

        sin(180 - x) = cos(x)
        cos(180 - x) = sin(x)
        tan(180 - x) = cot(x)
        cot(180 - x) = tan(x)
        sec(180 - x) = csc(x)
        csc(180 - x) = sec(x)

        sin(a + b) = sin(a) * cos(b) + sin(b) * cos(a)
        sin(a - b) = sin(a) * cos(b) - sin(b) * cos(a)
        cos(a + b) = cos(a) * cos(b) - sin(a) * sin(b)
        cos(a - b) = cos(a) * cos(b) + sin(a) * sin(b)
        tan(a + b) = (tan(a) + tab(b)) / (1 - tan(a) * tan(b))
        tan(a + b) = (tan(a) - tab(b)) / (1 + tan(a) * tan(b))

    Point and Vector:
        Vector: magnitude + angle
        Given (x, y) -> magnitude: sqrt(x^2 + y^2), angle: atan(y / x)
        Vector addition: (x1, y1) + (x2, y2) = (x1 + x2, y1 + y2)
        Vector subtraction: C = A - B -> A = B + C, (xc, yc) = (xa - xb, ya - yb)
        Dot Product:
            a.b = xa * xb + ya * yb + za * zb + ... = |a| * |b| * cos(theta)
            theta < 90: a.b > 0
            theta = 90: a.b = 0
            theta > 90: a.b < 0
            a.b = b.a
            a.(b + c) = a.b + a.c
            (c1 * a).(c2 * b) = c1 * c2 * (a.b)
            a.b = a.c -> a.(b - c) = 0 -> a, b-c are perpendicular
        Cross Product
            a * b = |a| * |b| * sin(theta)
            a, b are convex -> a * b > 0
            theta = 0 or 180 -> a * b = 0
            a, b are concave -> a * b < 0
            a * b = -b * a
            a * (b + c) = a * b + a * c
            a * b * c + b * c * a + c * a * b = 0
            |a * b|^2 = (|a| * |b|)^2 - (a.b)^2
        Standard Basis:
            ex = (1, 0), ey = (0, 1)
            i = ex = (1, 0, 0), j = ey = (0, 1, 0), k = ez = (0, 0, 1)
            i = j * k, j = k * i, k = i * j
            i * i = j * j = k * k = 0
            a * b = [i j k; xa ya za; xb yb zb]
        Euclidean Transformation:
            Translation:
                x' = x + h, y' = y + k
                [x'; y'; 1] = [1, 0, h; 0, 1, k; 0, 0, 1] * [x; y; 1]
            Rotation by angle T clockwise:
                x' = cos(T) * x - sin(T) * y
                y' = sin(T) * x + cos(T) * y
                [x'; y'; 1] = [cos(T), -sin(T), 0; sin(T), cos(T), 0; 0, 0, 1] * [x; y; 1]
            Rotation and Translation:
                [x'; y'; 1] = [cos(T), -sin(T), h; sin(T), cos(T), k; 0, 0, 1] * [x; y; 1]
            Translation and Rotation:
                [x'; y'; 1] = [cos(T), -sin(T), h*cos(T)-k*sin(T); sin(T), cos(T), h*sin(T)+k*cos(T); 0, 0, 1] * [x; y; 1]

    Complex Numbers and 2D Points
        Complex Numbers:
            Imaginary Unit i: i^2 = -1
            Complex Number: a + b*i (a, b are real)
        We can describe complex numbers in the coordinate plane as point (a, b)
        z = a + b*i = (a, b) = (r*cos(T), r*sin(T))
        r = sqrt(a^2 + b^2), T = atan(a / b)
        T = 0 -> z = 1
        T = pi/2 -> z = i
        T = pi -> z = -1
        T = 3pi/2 -> z = -i

        v1 = x1 + i*y1 -> conj(v1) = x1 - i*y1
        v2 = x2 + i*y2
        conj(v1) * (v2) = (x1*x2 + y1*y2) + (x1*y2 - x2*y1)*i = (dot product) + (cross product)*i

        revisit rotation and reflection

    Lines and Distances

*/

#include <bits/stdc++.h>
#include <complex>
#include <cmath>
#define int long long
#define double long double

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
#define square(a) (dot((a),(a)))

#define rotate0(a,ang) ((a)*exp(point(0,ang)))
#define rotate_about(a,ang,about) (rotate0(vec(about,a),ang)+about)

#define reflect0(a,b) (conj((a)/(b))*(b))
#define reflect_about(p, p0, p1) (reflect0(vec(p0, p), vec(p0, p1)) + p0)

double to_radians(double degree) { return degree / 180.0 * pi; }
double to_degree(double radian) { return radian / pi * 180.0; }
int dcmp(double x, double y) {
    if (fabs(x - y) < eps) return 0;
    return (x < y ? -1 : 1);
}

bool is_collinear(point a, point b, point c) {
    return dcmp(fabs(cross(b - a, c - a)), 0) == 0;
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

pair<double, point > circle(point a, point b, point c) {
    if (is_collinear(a, b, c))
        return {-1, point(0, 0)};

    point m1 = (a + b) / (double) 2.0, m2 = (a + c) / (double) 2.0;
    point v1 = b - a, v2 = c - a;
    point p1(-v1.Y, v1.X), p2(-v2.Y, v2.X);
    point end1 = m1 + p1, end2 = m2 + p2;
    point center = lines_intersect(m1, end1, m2, end2);
    return {length(a - center), center};
}
vector<point > circle_line_intersection(double r, point c, point a, point b) // radius, center, line ab
{
    double A = dot(b - a, b - a);
    double B = 2 * dot(b - a, a - c);
    double C = dot(a - c, a - c) - r * r;
    double disc = B * B - 4 * A * C;

    vector<point > res;
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
vector<point > circle_circle_intersection(double r1, point c1, double r2, point c2) {
    // same circle with positive radius
    if (same(c1, c2) && dcmp(r1, r1) == 0 && dcmp(r1, 0) == 1)
        return vector<point >(3, c1);

    auto get_angle = [&](double a, double b, double c) {
        return acos((b * b + c * c - a * a) / (2 * b * c));
    };

    double ang1 = angle(c2 - c1), ang2 = get_angle(r2, r1, length(c2 - c1));
    if (isnan(ang2))
        ang2 = 0;

    vector<point > res;
    point p1 = c1 + polar(r1, ang1 + ang2);

    if (dcmp(dot(p1 - c1, p1 - c1), r1 * r1) != 0 || dcmp(dot(p1 - c2, p1 - c2), r2 * r2) != 0)
        return vector<point >();

    res.push_back(p1);
    point p2 = c1 + polar(r1, ang1 - ang2);
    if (!same(p1, p2))
        res.push_back(p2);

    return res;
}

bool cmp(point a, point b) {
    double ang1 = to_degree(angle(a));
    double ang2 = to_degree(angle(b));

    if (dcmp(ang1, 0) == -1)
        ang1 += rotation;
    if (dcmp(ang2, 0) == -1)
        ang2 += rotation;

    if (dcmp(ang1, ang2) == 0)
        return false;
    return ang1 < ang2;
}

bool is_polygon_simple(vector<point > &p) {
    bool ok = true;
    p.push_back(p[0]);
    for (int i = 0; i < (int) p.size() - 1; i++) {
        for (int j = i + 2; j < (int) p.size() - 1; j++) {
            if (are_segments_intersect(p[i], p[i + 1], p[j], p[j + 1]) && !(i == 0 && j == (int) p.size() - 2))
                ok = false;
        }
    }
    p.pop_back();
    return ok;
}
bool is_polygon_convex(vector<point > &p) {
    bool ok = true;
    int sign = counter_clockwise(p[0], p[1], p[2]);
    p.push_back(p[0]);
    p.push_back(p[1]);

    for (int i = 1; i < (int) p.size() - 2; i++) {
        if (counter_clockwise(p[i], p[i + 1], p[i + 2]) != sign)
            ok = false;
    }

    p.pop_back();
    p.pop_back();
    return ok;
}
double polygon_area(vector<point > &p) {
    double area = 0;
    p.push_back(p[0]);
    for (int i = 0; i < (int) p.size() - 1; i++)
        area += cross(p[i], p[i + 1]);
    p.pop_back();
    return fabs(area / 2.0);
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

        if (are_segments_intersect(a, b, p[i], p[i + 1])) {
            point j = lines_intersect(a, b, p[i], p[i + 1]);
            right.push_back(j);
            left.push_back(j);
        }

        if (dcmp(cross(b - a, p[i] - a), 0) <= 0)
            left.push_back(p[i]);
    }
    return {left, right};
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
