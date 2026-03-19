#include <cstdio>
#include <fstream>
#define MAX_LENGTH 50

typedef struct msg {
    int p; // priority
    double tt; // transmission time;
    double T; // period
} Msg;

void print_info(int n, double tau, Msg m[]) {
    printf("%d, %f\n", n, tau);
    for (int i = 0; i < n; i++) {
        printf("%d %f %f\n", m[i].p, m[i].tt, m[i].T);
    }
}

double wcrt(int n, double tau, Msg m[]) { // worst case response time
    double res;
    // TODO
    return res;
}

int main() {
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

}