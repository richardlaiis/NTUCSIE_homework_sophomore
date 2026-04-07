sudo pacman -S intel-oneapi-mkl, openblas, cblas

gcc my_steepest_descent.c -O2 -Wall -Wextra -std=c11 -lmkl_rt -lm -o my_steepest_descent

gcc my_steepest_descent.c -O2 -Wall -Wextra -std=c11 -lopenblas -lm -o my_steepest_descent

sudo pacman -S intel-oneapi-mkl
gcc my_steepest_descent.c -O2 -Wall -Wextra -std=c11 -lcblas -lm -o my_steepest_descent