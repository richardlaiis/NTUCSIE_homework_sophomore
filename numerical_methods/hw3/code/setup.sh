sudo pacman -S intel-oneapi-mkl, openblas, cblas
source /opt/intel/oneapi/setvars.sh

gcc my_steepest_descent.c -O2 -Wall -Wextra -std=c11 -lmkl_rt -lm -o my_steepest_descent

gcc my_steepest_descent.c -O2 -Wall -Wextra -std=c11 -lopenblas -lm -o my_steepest_descent

sudo pacman -S intel-oneapi-mkl
gcc my_steepest_descent.c -O2 -Wall -Wextra -std=c11 -lcblas -lm -o my_steepest_descent

# https://github.com/amd/blis/blob/master/docs/BuildSystem.md
git clone https://github.com/flame/blis.git
cd blis/
./configure --enable-cblas auto
make
sudo make install
gcc -I /usr/local/include/blis my_steepest_descent.c -O2 -Wall -Wextra -std=c11 -o my_steepest_descent /usr/local/lib/libblis.a -lm
