#include <cstdio>
#include <fstream>
#include <cmath>
#include <chrono>
#define MAX_LENGTH 50
#define TIME_LIMIT 15000000 // in microseconds (1s = 1e6 microseconds)

typedef struct msg {
    int p; // priority
    double tt; // transmission time;
    double T; // period
} Msg;

double wcrt(int i, int n, double tau, Msg m[]);

void print_info(int n, double tau, Msg m[]) {
    printf("%d, %f\n", n, tau);
    for (int i = 0; i < n; i++) {
        printf("(%d %f %f): %f\n", m[i].p, m[i].tt, m[i].T, wcrt(i, n, tau, m));
    }
}

double wcrt(int i, int n, double tau, Msg m[]) { // worst case response time
    double Ci = m[i].tt;
    double Bi = 0.0;
    for (int k = 0; k < n; k++) {
        if (m[k].p >= m[i].p) {
            Bi = std::max(Bi, m[k].tt);
        }
    }
    
    double Qi = Bi;
    while (1) {
        double temp = 0.0;
        for (int k = 0; k < n; k++) {
            if (m[k].p < m[i].p) {
                temp += std::ceil((Qi+tau)/m[k].T)*m[k].tt;
            }
        }
        temp += Bi;

        if (temp > m[i].T) {
            return -1.0; // not schedulable
        }

        if (std::fabs(Qi-temp) < 1e12) {
            return Qi + Ci;
        } else {
            Qi = temp;
        }
    }
}

double sum_wcrt(int n, double tau, Msg m[]) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double Ri = wcrt(i, n, tau, m);
        if (Ri < 0) return -1.0;
        sum += Ri;
    }
    return sum;
}

int main() {
    auto beg = std::chrono::high_resolution_clock::now();

    int n;
    double tau;
    Msg m[MAX_LENGTH];
    int bestP[MAX_LENGTH]; // best priority assignment
    double minsum = __DBL_MAX__;

    std::ifstream fin("input.dat");
    fin >> n;
    fin >> tau;

    for (int i = 0; i < n; i++) {
        fin >> m[i].p >> m[i].tt >> m[i].T;
        bestP[i] = m[i].p;
    }

    print_info(n, tau, m);

    // simple timer
    while (1) {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
        if (duration.count() >= TIME_LIMIT) {
            printf("%d ms passed", duration.count());
            break;
        }
    }
    
}