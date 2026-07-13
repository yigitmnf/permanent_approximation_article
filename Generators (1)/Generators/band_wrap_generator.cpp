#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <n> <k> <output_path>\n";
    return 1;
  }

  const int n = std::atoi(argv[1]);
  const int k = std::atoi(argv[2]);
  const std::string output_path = argv[3];

  if (n <= 0) {
    std::cerr << "n must be > 0.\n";
    return 1;
  }
  if (k < 0) {
    std::cerr << "k must be >= 0.\n";
    return 1;
  }

  const int wrap_threshold = n - 1 - k;
  long long nnz = 0;
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= n; ++j) {
      if (std::abs(i - j) <= k || std::abs(i + j - (n + 1)) <= k) {
        nnz++;
      }
    }
  }

  std::ofstream out(output_path);
  if (!out) {
    std::cerr << "Failed to open output file: " << output_path << "\n";
    return 1;
  }

  out << "%%MatrixMarket matrix coordinate integer general\n";
  out << "% band/wrap matrix with ones for |i-j| <= k or |i-j| >= n-k\n";
  out << "% n=" << n << ", k=" << k << "\n";
  out << n << " " << n << " " << nnz << "\n";

  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= n; ++j) {
      if (std::abs(i - j) <= k || std::abs(i + j - (n + 1)) <= k) {
        out << i << " " << j << " 1\n";
      }
    }
  }

  return 0;
}
