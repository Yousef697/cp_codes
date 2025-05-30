#include <bits/stdc++.h>
#define int long long

using namespace std;

const double pi = acos(-1);
typedef complex<long double> point;

vector<point> fft(vector<point> P) {
    int n = P.size();
    if (n == 1) {
        return P;
    }

    vector<point> P_even(n / 2), P_odd(n / 2);
    for (int i = 0; i < n / 2; i++) {
        P_even[i] = P[2 * i];
        P_odd[i] = P[2 * i + 1];
    }
    vector<point> val_odd = fft(P_odd);
    vector<point> val_even = fft(P_even);
    vector<point> P_val(n);

    for (int i = 0; i < n / 2; i++) {
        point wj(cos(2 * pi * i / n), sin(2 * pi * i / n));
        P_val[i] = val_even[i] + wj * val_odd[i];
        P_val[i + n / 2] = val_even[i] - wj * val_odd[i];
    }

    return P_val;
}

vector<point> ifft(vector<point> P_val) {
    int n = P_val.size();
    if (n == 1) {
        return P_val;
    }

    vector<point> val_even(n / 2), val_odd(n / 2);
    for (int i = 0; i < n / 2; i++) {
        point wj(cos(2 * pi * i / n), sin(2 * pi * i / n));
        val_even[i] = (long double)0.5 * (P_val[i] + P_val[i + n / 2]);
        val_odd[i] = (long double)0.5 * (P_val[i] - P_val[i + n / 2]) / wj;
    }

    vector<point> P_even = ifft(val_even);
    vector<point> P_odd = ifft(val_odd);
    vector<point> P(n);
    for (int i = 0; i < n / 2; i++) {
        P[i * 2] = P_even[i];
        P[i * 2 + 1] = P_odd[i];
    }

    return P;
}

vector<int> multiply(vector<int> P, vector<int> Q) {
    int n = 1;
    while (n < P.size() + Q.size() - 1)
        n *= 2;

    vector<point> P_points(n), Q_points(n);
    copy(P.begin(), P.end(), P_points.begin());
    copy(Q.begin(), Q.end(), Q_points.begin());
    vector<point> P_val = fft(P_points);
    vector<point> Q_val = fft(Q_points);

    vector<int> R(n);
    vector<point> R_val(n);
    for (int i = 0; i < n; i++)
        R_val[i] = P_val[i] * Q_val[i];
    vector<point> R_points = ifft(R_val);

    for (int i = 0; i < n; i++) {
        R[i] = round(R_points[i].real());
    }

    return R;
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n, m;

    cin >> n;
    vector<int> a(n);
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }

    cin >> m;
    vector<int> b(m);
    for (int i = 0; i < m; i++) {
        cin >> b[i];
    }

    auto ans = multiply(a, b);

    for (int i = 0; i < n + m - 1; i++)
        cout << ans[i] << " ";

    return 0;
}
