#include <bits/stdc++.h>
using namespace std;

signed main() {
	ifstream fin("input.dat");
	int n;
	double tau;
	fin >> n;
	fin >> tau;
	vector<int> priority(n);
	vector<double> transmission_time(n), period(n);
	for (int i = 0; i < n; ++i) {
		fin >> priority[i] >> transmission_time[i] >> period[i];
	}
	fin.close();

    vector<double> res(n);
    for (int i = 0; i < n; i++) {
        double block = 0;
        for (int j = 0; j < n; j++) {
            if (priority[j] >= priority[i]) {
                block = max(block, transmission_time[j]);
            }
        }
        
        double wait = block;
        while (1) {
            double result  = block;
            for (int j = 0; j < n; j++) {
                if (priority[j] < priority[i]) {
                    result += transmission_time[j] * ceil((wait + tau) / period[j]); 
                }
            }
            if (result > period[i]) {
                res[i] = -1;
                break;
            }
            if (result-wait <= 1e-9 && result-wait >= -(1e-9)) {
                res[i] = transmission_time[i] + wait;
                break;
            } else {
                wait = result;
                continue;
            }
        }
    }

    for (int i = 0; i < n; i++) {
        cout << res[i] << endl;
    }
}