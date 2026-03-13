#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <ctime>

using namespace std;

using cd = complex<double>;
const double PI = acos(-1.0);


void fft(vector<cd>& a, bool invert) {
    int n = (int)a.size();

    
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j >= bit; bit >>= 1) j -= bit;
        j += bit;
        if (i < j) swap(a[i], a[j]);
    }

    
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        cd wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i + j];
                cd v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (cd& x : a) x /= n;
    }
}


vector<double> multiply(const vector<double>& a, const vector<double>& b) {
    int n1 = (int)a.size(), n2 = (int)b.size();
    int n = 1;
    while (n < n1 + n2 - 1) n <<= 1;          

    vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    fa.resize(n);
    fb.resize(n);

    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++) fa[i] *= fb[i];
    fft(fa, true);

    vector<double> result(n1 + n2 - 1);
    for (int i = 0; i < n1 + n2 - 1; i++) {
        result[i] = round(fa[i].real());      
    }
    return result;
}


vector<double> naive_multiply(const vector<double>& a, const vector<double>& b) {
    int n1 = (int)a.size(), n2 = (int)b.size();
    vector<double> result(n1 + n2 - 1, 0.0);
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            result[i + j] += a[i] * b[j];
        }
    }
    return result;
}

void print_poly(const vector<double>& p) {
    cout << "Coefficients (lowest degree → highest): ";
    for (double c : p) {
        cout << fixed << setprecision(0) << c << " ";
    }
    cout << endl;
}

int main() {
    cout << "=== FFT Polynomial Multiplication Demo ===\n\n";

    int d1, d2;
    cout << "Degree of first polynomial: ";
    cin >> d1;
    vector<double> p1(d1 + 1);
    cout << "Enter " << d1 + 1 << " coefficients (constant → x^" << d1 << "): ";
    for (auto& c : p1) cin >> c;

    cout << "Degree of second polynomial: ";
    cin >> d2;
    vector<double> p2(d2 + 1);
    cout << "Enter " << d2 + 1 << " coefficients: ";
    for (auto& c : p2) cin >> c;

    
    clock_t start = clock();
    vector<double> res_fft = multiply(p1, p2);
    clock_t end = clock();
    double time_fft = (end - start) / (double)CLOCKS_PER_SEC;

    start = clock();
    vector<double> res_naive = naive_multiply(p1, p2);
    end = clock();
    double time_naive = (end - start) / (double)CLOCKS_PER_SEC;

    cout << "\n=== FFT Result ===\n";
    print_poly(res_fft);
    cout << "Time: " << fixed << setprecision(6) << time_fft << " seconds\n";

    cout << "\n=== Naïve Result ===\n";
    print_poly(res_naive);
    cout << "Time: " << fixed << setprecision(6) << time_naive << " seconds\n";


    bool match = (res_fft.size() == res_naive.size());
    if (match) {
        for (size_t i = 0; i < res_fft.size(); ++i) {
            if (abs(res_fft[i] - res_naive[i]) > 1e-6) {
                match = false;
                break;
            }
        }
    }
    cout << "\nResults match: " << (match ? "YES" : "NO") << endl;

    return 0;
}