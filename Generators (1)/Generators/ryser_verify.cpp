#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

bool read_mtx(const std::string &path, int &n, std::vector<std::vector<int>> &A) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    bool symmetric = false;
    bool pattern = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '%') {
            if (line.find("symmetric") != std::string::npos) symmetric = true;
            if (line.find("pattern") != std::string::npos) pattern = true;
            continue;
        }
        std::istringstream iss(line);
        int n1, n2, nz;
        if (!(iss >> n1 >> n2 >> nz)) return false;
        n = n1;
        A.assign(n, std::vector<int>(n, 0));
        break;
    }
    int r, c, v;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '%') continue;
        std::istringstream iss(line);
        if (pattern) {
            if (!(iss >> r >> c)) continue;
            v = 1;
        } else {
            if (!(iss >> r >> c >> v)) continue;
        }
        if (r < 1 || c < 1 || r > n || c > n) continue;
        A[r-1][c-1] = v;
        if (symmetric && r != c) {
            A[c-1][r-1] = v;
        }
    }
    return true;
}

cpp_int ryser_permanent(const std::vector<std::vector<int>> &A) {
    int n = A.size();
    if (n == 0) return cpp_int(1);
    cpp_int total = 0;
    int max_mask = 1 << n;
    for (int mask = 1; mask < max_mask; ++mask) {
        int bits = __builtin_popcount(mask);
        cpp_int prod = 1;
        for (int i = 0; i < n; ++i) {
            int row_sum = 0;
            for (int j = 0; j < n; ++j) {
                if (mask & (1 << j)) row_sum += A[i][j];
            }
            prod *= row_sum;
            if (prod == 0) break;
        }
        if ((n - bits) % 2 == 1) total -= prod;
        else total += prod;
    }
    return total;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <matrix.mtx>\n";
        return 1;
    }
    int n;
    std::vector<std::vector<int>> A;
    if (!read_mtx(argv[1], n, A)) {
        std::cerr << "Failed to read " << argv[1] << "\n";
        return 1;
    }
    cpp_int perm = ryser_permanent(A);
    std::cout << perm << "\n";
    return 0;
}
