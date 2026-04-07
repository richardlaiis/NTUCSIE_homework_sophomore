sudo pacman -S intel-oneapi-mkl, openblas, cblas
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
source /opt/intel/oneapi/setvars.sh

cd ~
git clone https://github.com/flame/blis.git
cd blis/
./configure --enable-cblas auto
make
sudo make install
export BLIS_NUM_THREADS=1
cd ../code/

gcc -I /usr/local/include/blis my_steepest_descent.c -O2 -Wall -Wextra -std=gnu11 -o my_steepest_descent /usr/local/lib/libblis.a -lm

gcc my_steepest_descent.c -O2 -Wall -Wextra -std=gnu11 -lmkl_rt -lm -o my_steepest_descent

gcc my_steepest_descent.c -O2 -Wall -Wextra -std=gnu11 -lopenblas -lm -o my_steepest_descent

gcc my_steepest_descent.c -O2 -Wall -Wextra -std=gnu11 -lcblas -lm -o my_steepest_descent

curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
uv venv matrix
source matrix/bin/activate
uv pip install scipy, numpy