#include <stdio.h>
#include <stdint.h>
#include <string.h>

void print_bits(double d, const char* label) {
    uint64_t bits;
    memcpy(&bits, &d, sizeof(bits));
    printf("%s: %f | Raw Hex: 0x%016lx\n", label, d, bits);
}

int main() {
    double pos_zero = 0.0, neg_zero = -0.0;

    print_bits(pos_zero, "Literal 0.0   ");
    print_bits(neg_zero, "Literal -0.0  ");

    if (pos_zero == neg_zero) {
        printf("\nResult: The == operator says they are equal.\n");
    }
    
    printf("1 /  0.0 = %f\n", 1.0 / pos_zero);
    printf("1 / -0.0 = %f\n", 1.0 / neg_zero);

    return 0;
}