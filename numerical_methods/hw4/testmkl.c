#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
int main() {

}

/*
source /opt/intel/oneapi/setvars.sh

gcc -fopenmp -m64 -I"${MKLROOT}/include" testmkl.c \
    -L"${MKLROOT}/lib/intel64" \
    -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
    -lpthread -lm -ldl -o testmkl
*/