#include <stdio.h>
#include <float.h>

double max1(double x, double y) {
    return (x > y) ? x : y;
}

double max2(double x, double y) {
    return (x < y) ? y : x;
} 

int main() {
    double a = -5, b = -0.0;
    double c1 = a / max1(b, 0.0);
    double c2 = a / max2(b, 0.0);
    printf("%f %f\n", c1, c2);
}