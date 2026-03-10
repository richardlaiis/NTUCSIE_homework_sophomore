#include <stdio.h>
#include <math.h>

int main() {
    double values[] = {-0.0, 0.0, -1.0 / INFINITY, 1.0 / INFINITY, -1.0 * 0.0, 1.0 * 0.0, copysign(0.0, 1.0), copysign(0.0, -1.0)};

    const char *names[] = {"Literal -0.0", "Literal 0.0", "-1.0 / INFINITY", "1.0 / INFINITY", "-1.0 * 0.0", "1.0 * 0.0", "copysign(0, +)", "copysign(0, -)"};

    int n = sizeof(values) / sizeof(values[0]);

    printf("Method             | Value | signbit()\n");
    printf("-------------------|-------|----------\n");

    for (int i = 0; i < n; i++) {
        printf("%-18s | %4.1f  | %d\n", names[i], values[i], signbit(values[i]));
    }

    printf("\n--- Validation through Division ---\n");
    printf("1.0 / +0.0 = %f\n", 1.0 / values[1]);
    printf("1.0 / -0.0 = %f\n", 1.0 / values[0]);

    return 0;
}