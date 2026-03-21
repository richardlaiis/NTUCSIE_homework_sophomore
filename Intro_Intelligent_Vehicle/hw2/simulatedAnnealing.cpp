#include <cstdio>
#include <fstream>
#include <cmath>
#include <chrono>
#include <cfloat>
#include <ctime>
#define MAX_LENGTH 50
#define TIME_LIMIT 15000000 
const double INF = 1e100;

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

        if (std::fabs(Qi-temp) < 1e-9) {
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
        if (Ri < 0) return INF;
        sum += Ri;
    }
    return sum;
}

int main() {
    auto beg = std::chrono::high_resolution_clock::now();
    srand(time(NULL));

    int n;
    double tau;
    Msg m[MAX_LENGTH];
    int bestP[MAX_LENGTH]; // best priority assignment
    double minsum = DBL_MAX;

    std::ifstream fin("input.dat");
    fin >> n;
    fin >> tau;

    for (int i = 0; i < n; i++) {
        fin >> m[i].p >> m[i].tt >> m[i].T;
        bestP[i] = m[i].p;
    }

    double curr = sum_wcrt(n, tau, m);
    if (curr < minsum) {
        minsum = curr;
        for (int i = 0; i < n; i++) bestP[i] = m[i].p;
    }

    double temperature = 5.0;
    const double temp_min = 1e-6;
    const double alpha = 0.997;
    const int loops = 300;
    
    while (temperature > temp_min) {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
        if (duration.count() >= TIME_LIMIT) {
            printf("%d ms passed", duration.count());
            break;
        }

        for (int i = 0; i < loops; i++) {
            int a = rand() % n, b = rand() % n;
            while (a == b) b = rand() % n;
            std::swap(m[a].p, m[b].p);

            double cand = sum_wcrt(n, tau, m);
            double delta = cand - curr;

            bool accept = false;
            if (delta <= 0.0) accept = true;
            else {
                double prob = std::exp(-delta / temperature);
                double x = (double) rand() / (RAND_MAX + 1.0);
                if (x <= prob) accept = true;
            }

            if (accept) {
                curr = cand;
                if (cand < minsum) {
                    minsum = cand;
                    for (int i = 0; i < n; i++) bestP[i] = m[i].p;
                }
            } else {
                std::swap(m[a].p, m[b].p);
            }
        }

        temperature *= alpha;
    }

    for (int i = 0; i < n; i++) {
        printf("%d\n", bestP[i]);
    }
    printf("Optimized result: %f\n", minsum);
    
}