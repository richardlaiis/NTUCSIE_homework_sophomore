#include <bits/stdc++.h>
using namespace std;

double calc1(double f, double A, double B) {
    return 1.0 - 1.0 / (exp(A*f+B)+1.0);
}

double calc2(double f, double A, double B) {
    return exp(A*f+B) / (1.0+exp(A*f+B));
}

int main() {
    cout << calc1(1, -64, 0) << endl;
    // output: 0

    cout << calc2(1, -64, 0) << endl;
    // output: 1.60381e-28
}