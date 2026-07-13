#include <bits/stdc++.h>
#include <chrono>
#include <random>

#define EPSILON 0.000000001

using namespace std;
using namespace std::chrono;

template<typename T>
using Matrix = vector<vector<T>>;
int n; // matrix dimension
long double coef = 0;

long double *row_gammas;
long double *gamma_precomputed;
long double *dp;

long double *weight;
long double ub;

int NO_SAMPLES = 128 * 1024 * 1024;
std::string filename;
long double time_limit;
int period = 1024;
int md_bucket_scan_limit = 1;
bool use_fixed_seed = false;
unsigned int fixed_seed = 0;
unsigned int rng_stream_counter = 0;

chrono::steady_clock::time_point start_time;

std::mt19937 make_rng() {
	if(use_fixed_seed) {
		unsigned int stream_seed = fixed_seed + (1000003u * rng_stream_counter);
		rng_stream_counter++;
		return std::mt19937(stream_seed);
	}
	std::random_device rd;
	return std::mt19937(rd());
}

// Reads a Matrix Market (.mtx) file and returns the matrix in CSR format
void read_coords(const std::string &filename, Matrix<long double>& matrix, int& n) {   
    bool add_ji = false;
    bool pattern_only = false;
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);

    std::string line;
    // Skip comments
    while (std::getline(in, line)) {
        if(line.size() > 0 && line.rfind("%%MatrixMarket", 0) == 0 && line.find("symmetric") != std::string::npos) {
            add_ji = true;
        } 
        if(line.size() > 0 && line.rfind("%%MatrixMarket", 0) == 0 && line.find("pattern") != std::string::npos) {
            pattern_only = true;
        } 
        if (line.size() > 0 && line[0] != '%') break;
    }
    // std::cout << "Matrix is symmetric: " << add_ji << std::endl;
    // Read matrix size and number of nonzeros
    std::istringstream header(line);
    int M, N, NNZ;
    header >> M >> N >> NNZ;
    
	if(M != N) {
		std::cout << "Matrix is not square" << std::endl;
		exit(1);
	} else {
		n = M;
	}

	for (int i = 0; i < n; i++) {
		matrix.push_back(vector<long double>(n));
		for (int j = 0; j < n; j++) {
			matrix[i][j] = 0;
		}
	}

    int i, j;
    long double val;
    for (int k = 0; k < NNZ; ++k) {
        long double entry = 1.0L;
        if(pattern_only) {
            in >> i >> j;
        } else {
            in >> i >> j >> val;
            entry = val;
        }
		matrix[i-1][j-1] += entry;
		if(i != j && add_ji) {
            matrix[j-1][i-1] += entry;
        } 
    }
    in.close();
}

struct CSRMatrix {
    public:
        int n;
        vector<int> row_ptr;
        vector<int> col_idx;
        vector<long double> vals;
};

static string md_bucket_suffix() {
	if(md_bucket_scan_limit <= 1) return "";
	return "_B" + to_string(md_bucket_scan_limit);
}

static long double md_edge_survival_score(int current_vertex, int chosen_vertex,
										  const int* ptrs, const int* ids,
										  const int* degs, const bool* matched) {
	long double score = 1.0L;
	int delete_vertex = current_vertex;
	for(int rep = 0; rep < 2; rep++) {
		for(int ptr = ptrs[delete_vertex]; ptr < ptrs[delete_vertex + 1]; ptr++) {
			int nbr = ids[ptr];
			if(nbr == current_vertex || nbr == chosen_vertex || matched[nbr]) continue;
			int deg = degs[nbr];
			if(deg <= 1) return 0.0L;
			score *= ((long double)deg - 1.0L) / (long double)deg;
		}
		delete_vertex = chosen_vertex;
	}
	return score;
}

static long double md_combined_edge_weight(int vertex, int ptr, const CSRMatrix& A,
										   const CSRMatrix& At, int nnz) {
	if(vertex < n) return A.vals[ptr];
	return At.vals[ptr - nnz];
}

static void md_vertex_min_and_weighted_ess(int vertex, const CSRMatrix& A,
										   const CSRMatrix& At, int nnz,
										   const int* ptrs, const int* ids,
										   const int* degs, const bool* matched,
										   long double& min_score, long double& ess_score) {
	bool has_option = false;
	min_score = numeric_limits<long double>::infinity();
	long double total_weight = 0.0L;
	long double sum = 0.0L;
	long double sumsq = 0.0L;
	for(int ptr = ptrs[vertex]; ptr < ptrs[vertex + 1]; ptr++) {
		int nbr = ids[ptr];
		if(matched[nbr]) continue;
		long double edge_weight = md_combined_edge_weight(vertex, ptr, A, At, nnz);
		total_weight += edge_weight;
		has_option = true;
	}
	if(!has_option || total_weight <= 0.0L) {
		min_score = 0.0L;
		ess_score = 0.0L;
		return;
	}
	for(int ptr = ptrs[vertex]; ptr < ptrs[vertex + 1]; ptr++) {
		int nbr = ids[ptr];
		if(matched[nbr]) continue;
		long double survival = md_edge_survival_score(vertex, nbr, ptrs, ids, degs, matched);
		long double edge_weight = md_combined_edge_weight(vertex, ptr, A, At, nnz);
		long double option_score = (edge_weight / total_weight) * survival;
		if(option_score < min_score) min_score = option_score;
		sum += option_score;
		sumsq += option_score * option_score;
	}
	ess_score = (sum > 0.0L && sumsq > 0.0L) ? (sum * sum) / sumsq : 0.0L;
}

static int md_take_vertex_from_bucket_weighted(int degree, int* q_ptrs_start, int* q_ptrs_end,
											  int* que, int* que_locs, int* degs,
											  bool* matched, int* ptrs, int* ids,
											  const CSRMatrix& A, const CSRMatrix& At,
											  int nnz) {
	while(q_ptrs_start[degree] < q_ptrs_end[degree] && degs[que[q_ptrs_start[degree]]] == 0) {
		q_ptrs_start[degree]++;
	}
	if(q_ptrs_start[degree] == q_ptrs_end[degree]) return -1;
	if(degree <= 1 || md_bucket_scan_limit <= 1) return que[q_ptrs_start[degree]++];

	int best_loc = -1;
	long double best_min_score = -1.0L;
	long double best_ess_score = -1.0L;
	int scanned = 0;
	for(int loc = q_ptrs_start[degree]; loc < q_ptrs_end[degree] && scanned < md_bucket_scan_limit; loc++) {
		int vertex = que[loc];
		if(degs[vertex] != degree || matched[vertex]) continue;
		scanned++;
		long double min_score = 0.0L;
		long double ess_score = 0.0L;
		md_vertex_min_and_weighted_ess(vertex, A, At, nnz, ptrs, ids, degs, matched, min_score, ess_score);
		if(min_score > best_min_score ||
		   (min_score == best_min_score && ess_score > best_ess_score) ||
		   (min_score == best_min_score && ess_score == best_ess_score && (best_loc < 0 || vertex < que[best_loc]))) {
			best_min_score = min_score;
			best_ess_score = ess_score;
			best_loc = loc;
		}
	}

	if(best_loc < 0) return que[q_ptrs_start[degree]++];

	int head_loc = q_ptrs_start[degree];
	int current_vertex = que[best_loc];
	if(best_loc != head_loc) {
		int head_vertex = que[head_loc];
		que[head_loc] = current_vertex;
		que[best_loc] = head_vertex;
		que_locs[current_vertex] = head_loc;
		que_locs[head_vertex] = best_loc;
	}
	q_ptrs_start[degree]++;
	return current_vertex;
}

CSRMatrix build_transpose(const CSRMatrix &A) {
    int n = A.n;
    std::vector<std::vector<std::pair<int, long double>>> adj(n);

    for (int i = 0; i < n; ++i) {
        for (int k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            int j = A.col_idx[k];
            adj[j].push_back({i, A.vals[k]});
        }
    }
    
    CSRMatrix B;
    B.n = n;
    B.row_ptr.resize(n + 1, 0);
    int nnz = 0;
    for (int i = 0; i < n; ++i) {
        B.row_ptr[i + 1] = B.row_ptr[i] + adj[i].size();
        nnz += adj[i].size();
    }
    B.col_idx.reserve(nnz);
    B.vals.reserve(nnz);
    for (const auto &lst : adj) {
        for (const auto &edge : lst) {
            B.col_idx.push_back(edge.first);
            B.vals.push_back(edge.second);
        }
    }
    return B;
}

void convertDenseToSparse(Matrix<long double>& matrix, CSRMatrix& A, int n) {
	A.n = n;

	int totalNNZ = 0;
	for (int i = 0; i < n; i++) {
		A.row_ptr.push_back(totalNNZ);
		for (int j = 0; j < n; j++) {
			if(matrix[i][j] > EPSILON) {
				A.col_idx.push_back(j);
				A.vals.push_back(matrix[i][j]);
				totalNNZ++;
			}
		}
	}
	A.row_ptr.push_back(totalNNZ);
}

template<typename T>
static void print_matrix_stats(const Matrix<T> &matrix, const string &label) {
	int n = (int)matrix.size();
	long long nnz = 0;
	long long max_nnz = 0;
	for (int i = 0; i < n; i++) {
		long long row_nnz = 0;
		for (int j = 0; j < n; j++) {
			if (matrix[i][j] > EPSILON) {
				nnz++;
				row_nnz++;
			}
		}
		if (row_nnz > max_nnz) max_nnz = row_nnz;
	}
	long double avg_nnz = n ? (long double)nnz / (long double)n : 0;
	cerr<<label<<": n="<<n<<" nnz="<<nnz<<" max_nnz="<<max_nnz<<" avg_nnz="<<fixed<<setprecision(2)<<avg_nnz<<endl;
}

template <typename T>
vector<int> hopcroft_karp(Matrix<T> matrix) {
	int n = matrix.size();
	vector<vector<int>> bipartite_edges(n);
	vector<int> matching(n, -1);
	vector<int> inv(n, -1);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (matrix[i][j] >= EPSILON) {
				bipartite_edges[i].push_back(j + n);
			}
		}
	}

	vector<vector<int>> dag(2 * n);
	bool modified = false;
	do {
		modified = false;
		queue<int> q;
		vector<int> seen(2 * n, 0);
		for (int i = 0; i < n; i++) {
			dag[i].clear();
			dag[i + n].clear();
			if (matching[i] == -1) {
				q.push(i);
			}
		}
		bool stop = false;
		while (!q.empty()) {
			int i = q.front();
			q.pop();
			if (seen[i]) continue;
			seen[i] = 1;
			if (i >= n) {
				if (inv[i - n] == -1) {
					stop = true;
				} else if (!seen[inv[i - n]] && !stop) {
					q.push(inv[i - n]);
					dag[i].push_back(inv[i - n]);
				}
			} else {
				for (int j : bipartite_edges[i]) {
					if (!seen[j] && !stop && matching[i] != j) {
						q.push(j);
						dag[i].push_back(j);
					}
				}
			}
		}

		vector<int> parent(2 * n, -1);
		vector<bool> handled(2 * n, false);
		stack<int> s;
		for (int i = 0; i < n; i++) {
			if (matching[i] == -1) {
				s.push(i);
				parent[i] = -2;
			}
			while (!s.empty()) {
				int i = s.top();
				s.pop();
				if (handled[i]) continue;
				handled[i] = true;
				if (i >= n && inv[i - n] == -1) {
					while (!s.empty()) {
						s.pop();
					}
					modified = true;
					int j = i;
					while (j != -2) {
						matching[parent[j]] = j;
						inv[j - n] = parent[j];
						j = parent[parent[j]];
					}
					break;
				} else {
					for (int j : dag[i]) {
						if (handled[j]) continue;
						parent[j] = i;
						s.push(j);
					}
				}
			}
		}
	} while (modified);

	vector<int> perm(n);
	for (int i = 0; i < n; i++) {
		perm[i] = matching[i] - n;
	}
	return perm;
}

template<typename T>
vector<int> scc(Matrix<T> m) {
	int n = m.size();
	vector<vector<int>> v(n);
	vector<vector<int>> r(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (m[i][j] > EPSILON) {
				v[i].push_back(j);
				r[j].push_back(i);
			}
		}
	}
	vector<int> state(n);
	vector<int> order;

	for (int i = 0; i < n; i++) {
		if (state[i]) continue;
		stack<int> s;
		s.push(i);
		state[i] = 1;
		while (!s.empty()) {
			int i = s.top();
			if (!v[i].empty()) {
				int j = v[i].back();
				v[i].pop_back();
				if (!state[j]) {
					s.push(j);
					state[j] = 1;
				}
			} else {
				s.pop();
				order.push_back(i);
			}
		}
	}
	vector<int> components(n);
	int component_id = 1;
	reverse(order.begin(), order.end());
	for (int i : order) {
		if (components[i]) continue;
		stack<int> s;
		s.push(i);
		while (!s.empty()) {
			int i = s.top();
			s.pop();
			if (components[i]) {
				continue;
			}
			components[i] = component_id;
			for (int j : r[i]) {
				s.push(j);
			}
		}
		component_id += 1;
	}

	return components;
}

template<typename T>
Matrix<T> tassa(Matrix<T> matrix) {
	int n = matrix.size();
	vector<int> matching = hopcroft_karp(matrix);
	int p[n];
	int pr[n];
	for (int i = 0; i < n; i++) {
		if (matching[i] < 0) return {};
		p[i] = matching[i];
		pr[matching[i]] = i;
		matching[i] = i;
	}
	{
		Matrix<T> permuted;
		for (int i = 0; i < n; i++) {
			permuted.push_back(matrix[pr[i]]);
		}
		matrix = permuted;
	}
	vector<long double> cs(n);
	for (int i = 0; i < n; i++) {
		cs[i] = matrix[i][i];
		matrix[i][i] = 0;
	}
	vector<int> components = scc(matrix);
	Matrix<T> support = matrix;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			support[i][j] = 0;
		}
	}
	for (int i = 0; i < n; i++) {
		support[i][i] = cs[i];
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (matrix[i][j] > EPSILON && components[i] == components[j]) {
				support[i][j] = matrix[i][j];
			}
		}
	}
	{
		Matrix<T> permuted;
		for (int i = 0; i < n; i++) {
			permuted.push_back(support[p[i]]);
		}
		support = permuted;
	}
	return support;
}

struct dataPoint {
	public:
		dataPoint(int _noSamples, long double _estimation, double _seconds) {
			noSamples = _noSamples;
			estimation = _estimation;
			seconds = _seconds;
		}

		int noSamples;
		long double estimation;
		double seconds;
}; 

template<typename ST>
void light_scale(CSRMatrix& A, CSRMatrix& At, ST* rv, ST* cv, bool* r_is_matched, bool* c_is_matched,  int current_vertex) {
	int n = A.n;

	if(current_vertex < n) { //this is a row vertex 
		ST row_sum = 0;
		for(int ptr = A.row_ptr[current_vertex]; ptr < A.row_ptr[current_vertex + 1]; ptr++) {
			int c = A.col_idx[ptr];
			if(c_is_matched[c] == false) {
				ST edge = static_cast<ST>(A.vals[ptr]);
				ST col_sum = 0;
				for(int ptr2 = At.row_ptr[c]; ptr2 < At.row_ptr[c + 1]; ptr2++) {
					int r = At.col_idx[ptr2];
					if(r_is_matched[r] == false) {

						// third level: nbr-nbr-row update using its neighbor columns
						ST row_sum2 = 0;
						for(int ptr3 = A.row_ptr[r]; ptr3 < A.row_ptr[r + 1]; ptr3++) {
							int c2 = A.col_idx[ptr3];
							if(c_is_matched[c2] == false) {
								ST edge2 = static_cast<ST>(A.vals[ptr3]);
								row_sum2 += edge2 * cv[c2];
							}
						}
						if(row_sum2 != 0) {
							rv[r] = 1/row_sum2;
						}

						col_sum += static_cast<ST>(At.vals[ptr2]) * rv[r];
					}
				}
				if(col_sum != 0) {
					cv[c] = 1/col_sum;
				}
				row_sum += edge * cv[c];
			}
		}
		if(row_sum != 0) {
			rv[current_vertex] = 1/row_sum;
		}
	} else {
		current_vertex -= n;
		ST col_sum = 0;
		for(int ptr = At.row_ptr[current_vertex]; ptr < At.row_ptr[current_vertex + 1]; ptr++) {
			int r = At.col_idx[ptr];
			if(r_is_matched[r] == false) {
				ST edge = static_cast<ST>(At.vals[ptr]);
				ST row_sum = 0;
				for(int ptr2 = A.row_ptr[r]; ptr2 < A.row_ptr[r + 1]; ptr2++) {
					int c = A.col_idx[ptr2];
					if(c_is_matched[c] == false) {

						// third level: nbr-nbr-col update using its neighbor rows
						ST col_sum2 = 0;
						for(int ptr3 = At.row_ptr[c]; ptr3 < At.row_ptr[c + 1]; ptr3++) {
							int r2 = At.col_idx[ptr3];
							if(r_is_matched[r2] == false) {
								ST edge2 = static_cast<ST>(At.vals[ptr3]);
								col_sum2 += edge2 * rv[r2];
							}
						}
						if(col_sum2 != 0) {
							cv[c] = 1/col_sum2;
						}

						row_sum += static_cast<ST>(A.vals[ptr2]) * cv[c];
					}
				}
				if(row_sum != 0) {
					rv[r] = 1/row_sum;
				}
				col_sum += edge * rv[r];
			}
		}
		if(col_sum != 0) {
			cv[current_vertex] = 1/col_sum;
		}
	}
}


template<typename ST, bool row_matching_given, bool col_matching_given, bool with_error, bool with_init>
void scale(	CSRMatrix& A, CSRMatrix& At, 
			ST* rv, ST* cv, 
			bool* r_is_matched, bool* c_is_matched, 
			int start_row, int start_col, int iter) {
	int n = A.n;

	if constexpr (with_init) {
		for(int i = 0; i < n; i++) {
			rv[i] = cv[i] = 1;
		}
	}

	for(int x = 0; x < iter; x++) {	
		for(int c = start_col; c < n; c++) {
			if constexpr (col_matching_given) {
				if(c_is_matched[c] == true) {
					continue;
				}
			}
			ST sum = 0;
			for(int ptr = At.row_ptr[c]; ptr < At.row_ptr[c + 1]; ptr++) {
				int r = At.col_idx[ptr];
				
				if(r >= start_row) {
					if constexpr (row_matching_given) {
						if(!r_is_matched[r]) {
							sum += static_cast<ST>(At.vals[ptr]) * rv[r];
						}
					} else {
						sum += static_cast<ST>(At.vals[ptr]) * rv[r];
					}
				}
			}
			if(sum != 0) {
				cv[c] = 1/sum;
			}
		}

			for(int r = start_row; r < n; r++) {
				if(r_is_matched == nullptr || r_is_matched[r] == false) {
					ST sum = 0;
					for(int ptr = A.row_ptr[r]; ptr < A.row_ptr[r + 1]; ptr++) {
						int c = A.col_idx[ptr];

						if(c >= start_col) {
							if constexpr (col_matching_given) {
								if(!c_is_matched[c]) {
									sum += static_cast<ST>(A.vals[ptr]) * cv[c];
								}
							} else {
								sum += static_cast<ST>(A.vals[ptr]) * cv[c];
							}
						}
					}
					if(sum != 0) {
					rv[r] = 1/sum;
				}
			}
		}
	}

	if constexpr (with_error) {
		long double max_error = 0;
		for(int c = start_col; c < n; c++) {
			if(c_is_matched == nullptr || c_is_matched[c] == false) {
				long double sum = 0;
				for(int ptr = At.row_ptr[c]; ptr < At.row_ptr[c + 1]; ptr++) {
					int r = At.col_idx[ptr];
					if(r >= start_row) {
						if constexpr (row_matching_given) {
							if(!r_is_matched[r]) {
								sum += static_cast<long double>(At.vals[ptr]) * cv[c] * rv[r];
							}
						} else {
							sum += static_cast<long double>(At.vals[ptr]) * cv[c] * rv[r];
						}
					}
				}
				if(sum != 0) {
					if(fabs(sum - 1) > max_error) {
						max_error = fabs(sum - 1);
					}
				}
			}
		}
		cout << "Max column error is: " << max_error << endl;

		max_error = 0;
		for(int r = start_row; r < n; r++) {
			if(r_is_matched == nullptr || r_is_matched[r] == false) {
				long double sum = 0;
				for(int ptr = A.row_ptr[r]; ptr < A.row_ptr[r + 1]; ptr++) {
					int c = A.col_idx[ptr];

					if(c >= start_col) {
						if constexpr (col_matching_given) {
							if(!c_is_matched[c]) {
								sum += static_cast<long double>(A.vals[ptr]) * cv[c] * rv[r];
							}
						} else {
							sum += static_cast<long double>(A.vals[ptr]) * cv[c] * rv[r];
						}
					}
				}
				if(sum != 0) {
					if(fabs(sum - 1) > max_error) {
						max_error = fabs(sum - 1);
					}
				}
			}
		}
		cout << "Max row error is: " << max_error << endl;
	}
}

//scale_approach is:
//0: no scaling
//1: once at the beginning
//2: once + scaling at every step
//3: once + light scaling at every step
//4: once + hybrid scaling at every step
//5: once + scaling at every step with interpolation
//6: once + scaling at every step with reverse interpolation
template<typename PT, typename RT, typename ST, int scale_approach, bool amplify>
void rasmussen_generic(CSRMatrix& A, CSRMatrix& At,
					   int no_outer_scaling_iter, int no_inner_scaling_iter, 
					   int hybrid_full_scale_period) {
	std::mt19937 gen = make_rng();
	std::uniform_real_distribution<RT> dis(0, 1);

	int n = A.n;
	vector<dataPoint> dataPoints; dataPoints.reserve(20000000);
	vector<PT> products; products.reserve(20000000);
	unsigned long* stats = new unsigned long[n + 1]; memset(stats, 0, sizeof(unsigned long) * (n + 1));

	bool *c_is_matched = new bool[n];
	bool *r_is_matched = new bool[n];
	int* unmatched_c = new int[n];
	RT* unmatched_terms = new RT[n];
	PT estimation = 0;

	ST *rv = nullptr, *cv = nullptr, *rv_sample = nullptr, *cv_sample = nullptr;
	if constexpr(scale_approach > 0) { //this means we do scaling at some point
		rv = new ST[n];
		cv = new ST[n];
		scale<ST, false, false, false, true>(A, At, rv, cv, nullptr, nullptr, 0, 0, no_outer_scaling_iter);
		
		if constexpr(scale_approach > 1) { //this means we do scaling at some point
			rv_sample = new ST[n];
			cv_sample = new ST[n];
		} else {
			rv_sample = rv;
			cv_sample = cv;
		}
	}

	auto start_timex = high_resolution_clock::now();
	int x = 0; //sample ID
	while (true) {//get a sample
		x++;
		memset(c_is_matched, false, sizeof(bool) * n); //all columns in the row that are unmatched
		memset(r_is_matched, false, sizeof(bool) * n); //all rows in the column that are unmatched

		if constexpr(scale_approach > 1) {
			memcpy(rv_sample, rv, sizeof(ST) * n);
			memcpy(cv_sample, cv, sizeof(ST) * n);
		}

		PT prod = 1; //estimation for sample - rv value
		int i = 0;
		for(; i < n; i++) { //for every row

			//inner scale logic
			if constexpr(scale_approach == 2) {
				scale<ST, false, true, false, false>(A, At, rv_sample, cv_sample, nullptr, c_is_matched, i, 0, no_inner_scaling_iter);
			} else if constexpr(scale_approach == 3) {
				light_scale<ST>(A, At, rv_sample, cv_sample, r_is_matched, c_is_matched, i);
			} else if constexpr(scale_approach == 4) {
				if(i % hybrid_full_scale_period == 0) {
					scale<ST, false, true, false, false>(A, At, rv_sample, cv_sample, nullptr, c_is_matched, i, 0, no_inner_scaling_iter);
				} else {
					light_scale<ST>(A, At, rv_sample, cv_sample, r_is_matched, c_is_matched, i);
				}
			} else if constexpr(scale_approach == 5) {
				scale<ST, false, true, false, false>(A, At, rv_sample, cv_sample, nullptr, c_is_matched, i, 0, (2 * (double(i)/n) * no_inner_scaling_iter) + 1);
			}  else if constexpr(scale_approach == 6) {
				scale<ST, false, true, false, false>(A, At, rv_sample, cv_sample, nullptr, c_is_matched, i, 0, (2 * no_inner_scaling_iter) - (2 * (double(i)/n) * no_inner_scaling_iter) + 1);
			}

			int no_unmatched = 0; //find number of unmatched columns
			RT sum_unmatched = 0;
			for(int ptr = A.row_ptr[i]; ptr < A.row_ptr[i+1]; ptr++) {
				int col_id = A.col_idx[ptr];
				if(c_is_matched[col_id] == false) {
					RT edge = static_cast<RT>(A.vals[ptr]);
					RT proposal_term = edge;
					if constexpr(scale_approach > 0) {
						if constexpr(amplify) {
							proposal_term *= cv_sample[col_id] * cv_sample[col_id];
						} else {
							proposal_term *= cv_sample[col_id];
						}
					}
					sum_unmatched += proposal_term;
					unmatched_terms[no_unmatched] = proposal_term;
					unmatched_c[no_unmatched++] = col_id;
				}
			}
			if(no_unmatched == 0) {
				prod = 0;
				break; //if everything is matched sample contribution is 0
			}

			int chosen_column = -1;
			RT rand_val = sum_unmatched * dis(gen);
			chosen_column = unmatched_c[0]; 
			RT running_sum = 0;
			for(int chosen = 0; chosen < no_unmatched; chosen++) {
				chosen_column = unmatched_c[chosen];
				running_sum += unmatched_terms[chosen];
				if(rand_val <= running_sum) {
					break;
				}
			}
			if constexpr(scale_approach == 0) {
				prod = prod * sum_unmatched;
			} else if constexpr(amplify) {
				prod = prod * (sum_unmatched / (cv_sample[chosen_column] * cv_sample[chosen_column]));
			} else {	
				prod = prod * (sum_unmatched / cv_sample[chosen_column]);
			}

			r_is_matched[i] = true;
			c_is_matched[chosen_column] = true;
		}

		stats[i]++;
		estimation += prod;
		products.push_back(prod);

		if(x % period == 0) {
			long double currentEstimation = estimation / x;
			auto duration = duration_cast<milliseconds>(high_resolution_clock::now() - start_timex);
			double seconds = ((double)duration.count()) / 1000;
			dataPoint p(x, currentEstimation, seconds);
			dataPoints.push_back(p);
			if(seconds > time_limit) {
				break;
			}
		}
	}
	auto stop_timex = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop_timex - start_timex);

	long int NO_SAMPLES = products.size();
	long double sum = 0;
	for(int i = 0; i < NO_SAMPLES; i++) {
		sum += products[i];
	}
	long double mean = sum / NO_SAMPLES;

	long double variance = 0;
	for(int i = 0; i < NO_SAMPLES; i++) {
		variance += (products[i] - mean) * (products[i] - mean);
	}
	variance /= (NO_SAMPLES - 1);
	long double stdev = sqrt(variance);

	long double stat_mean = 0;
	int q4pt = -1, q3pt = -1, q2pt = -1;
	long double sum_stat = 0;
	int max_cardinality = 0;
	for(long long i = 0; i <= n; i++) {
		sum_stat += stats[i];
		if((q4pt == -1) && ((sum_stat / NO_SAMPLES) > 0.25)) q4pt = i;
		if((q3pt == -1) && ((sum_stat / NO_SAMPLES) > 0.50)) q3pt = i;
		if((q2pt == -1) && ((sum_stat / NO_SAMPLES) > 0.75)) q2pt = i;
		if(stats[i] > 0) {
			max_cardinality = i;
		}
		stat_mean += i * stats[i];
	}
	stat_mean /= NO_SAMPLES;

	cout << "rasmussen_"; 
	if constexpr(scale_approach == 0) {
		cout << "no_scaling";
	} else if constexpr(scale_approach == 1) {
		cout << "once_scaling";
	} else if constexpr(scale_approach == 2) {
		cout << "scaling_every_step";
	} else if constexpr(scale_approach == 3) {
		cout << "light_scaling_every_step";
	} else if constexpr(scale_approach == 4) {
		cout << "hybrid_scaling_every_step";
	} else if constexpr(scale_approach == 5) {
		cout << "scaling_every_step_with_interpolation";
	} else if constexpr(scale_approach == 6) {
		cout << "scaling_every_step_with_rev_interpolation";
	}
	if(amplify) {
		cout << "_with_amplification";
	} 
	if(md_bucket_scan_limit > 1) {
		cout << "_B" << md_bucket_scan_limit;
	}
	
	cout << " with " << NO_SAMPLES << " samples ------------------------------------------" << endl;
	cout << "Estimation: " << mean << "\tSamples_variance: " << variance << "\tSamples_stdev: " << stdev << endl;
	cout << "Ratio_perfect_matching: " << ((long double)stats[n]) / NO_SAMPLES << "\tMean_matching_cardinality: " << stat_mean << "\tQs: [" 
		 << q4pt << ", " << q3pt << ", " << q2pt << "] -- maxcard:" << max_cardinality << endl;
    cout << "Total_time_sec: " << ((double)duration.count()) / 1000 << "\tTime_taken_by_a_single_sample_ms: " << ((double)duration.count()) / NO_SAMPLES << endl << endl;

	ofstream out_file;
	string run_bucket_suffix = md_bucket_suffix();
	if constexpr(scale_approach == 0) {
		if constexpr(amplify) {
			out_file.open (filename + ".noscl" + run_bucket_suffix + "_amplify.wout");
		} else {
	  		out_file.open (filename + ".noscl" + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 1) {
		if constexpr(amplify) {
			out_file.open (filename + ".oncescl" + run_bucket_suffix + "_amplify.wout");
		} else {
	  		out_file.open (filename + ".oncescl" + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 2) {
		if constexpr(amplify) {
			out_file.open (filename + ".allscl" + run_bucket_suffix + "_amplify.wout");
		} else {
  			out_file.open (filename + ".allscl" + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 3) {
		if constexpr(amplify) {
			out_file.open (filename + ".lightscl" + run_bucket_suffix + "_amplify.wout");
		} else {
  			out_file.open (filename + ".lightscl" + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 4) {
		if constexpr(amplify) {
			out_file.open (filename + ".hybridscl_P" + to_string(hybrid_full_scale_period) + run_bucket_suffix + "_amplify.wout");
		} else {
  			out_file.open (filename + ".hybridscl_P" + to_string(hybrid_full_scale_period) + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 5) {
		if constexpr(amplify) {
			out_file.open (filename + ".intscl_P" + to_string(hybrid_full_scale_period) + run_bucket_suffix + "_amplify.wout");
		} else {
  			out_file.open (filename + ".intscl_P" + to_string(hybrid_full_scale_period) + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 6) {
		if constexpr(amplify) {
			out_file.open (filename + ".revintscl_P" + to_string(hybrid_full_scale_period) + run_bucket_suffix + "_amplify.wout");
		} else {
  			out_file.open (filename + ".revintscl_P" + to_string(hybrid_full_scale_period) + run_bucket_suffix + ".wout");
		}
	}			
	out_file << n << "\n";
	for(int i = 1; i < n; i++) {
		out_file << stats[i] << " ";
	}
	out_file << stats[n] << "\n";
	out_file << dataPoints.size() << "\n";
	for(auto dp : dataPoints) {
		out_file << dp.noSamples << "\t" << dp.estimation << "\t" << dp.seconds << "\n";
	}
	out_file.close();

	delete [] c_is_matched;
	delete [] r_is_matched;

	delete [] unmatched_c;
	delete [] unmatched_terms;
	delete [] stats;
	if constexpr(scale_approach > 0) {
		delete [] rv;
		delete [] cv;
		if constexpr(scale_approach > 1) {
			delete [] rv_sample;
			delete [] cv_sample;
		}
	}
}

//scale_approach is:
//0: no scaling
//1: once at the beginning
//2: once + scaling at every step
//3: once + light scaling at every step
//4: once + hybrid scaling at every step
//5: once + scaling at every step with reverse interpolation
//7: once + adaptive scaling at every step
template<typename PT, typename RT, typename ST, int scale_approach, bool amplify, bool defer_degree_2 = false>
void rasmussen_md_generic(CSRMatrix& A, CSRMatrix& At,
								int no_outer_scaling_iter, int no_inner_scaling_iter, 
								int hybrid_full_scale_period) {
	std::mt19937 gen = make_rng();
	std::uniform_real_distribution<RT> dis(0, 1);

	int nnz = A.row_ptr[n];
	int* ptrs = new int[(2 * n) + 1];
	int* ids = new int[2 * nnz];

	for(int i = 0; i < n; i++) {
		ptrs[i] = A.row_ptr[i];
		ptrs[i + n] = At.row_ptr[i] + nnz;
	}
	ptrs[2 * n] = 2 * nnz;

	for(int i = 0; i < nnz; i++) {
		ids[i] = A.col_idx[i] + n;
		ids[i + nnz] = At.col_idx[i];
	}

	vector<dataPoint> dataPoints; dataPoints.reserve(100000000);
	vector<PT> products; products.reserve(100000000);

	unsigned long* stats = new unsigned long[n + 1];
	memset(stats, 0, sizeof(unsigned long) * (n + 1));

	int* degs = new int[2 * n];
	int* deg_counts = new int[n + 1];
	memset(deg_counts, 0, sizeof(int) * (n + 1));
	for(int i = 0; i < 2 * n; i++) {
		degs[i] = ptrs[i + 1] - ptrs[i];
		deg_counts[degs[i]]++;
	}

	int* q_ptrs_start = new int[n + 1];
	int* q_ptrs_end = new int[n + 1];
	q_ptrs_start[0] = 0;
	for(int i = 1; i <= n; i++) q_ptrs_start[i] = q_ptrs_start[i - 1] + deg_counts[i - 1];
	for(int i = 0; i <= n; i++) q_ptrs_end[i] = q_ptrs_start[i] + deg_counts[i];

	int* que = new int[2 * n];
	int *que_locs = new int[2 * n];
	for(int i = 0; i < 2 * n; i++) {
		int &loc = q_ptrs_start[degs[i]];
		que_locs[i] = loc;
		que[loc++] = i;
	}
	for(int i = n; i > 0; i--) q_ptrs_start[i] = q_ptrs_start[i-1];
	q_ptrs_start[0] = 0;

	int min_deg = 0;
	for(int i = 0; i <= n; i++) {
		if(deg_counts[i] > 0) {
			min_deg = i;
			break;
		} 
	}

	if(min_deg == 0) {
		cout << "Permanent is 0 - there is an empty row/column" << endl;
		exit(1); 
	}

	ST *rcv = nullptr, *rcv_sample = nullptr;
	if constexpr(scale_approach > 0) { //this means we do scaling at some point
		rcv = new ST[2 * n];
		scale<ST, false, false, false, true>(A, At, rcv, rcv + n, nullptr, nullptr, 0, 0, no_outer_scaling_iter);
		
		if constexpr(scale_approach > 1) { //this means we do scaling at some point
			rcv_sample = new ST[2 * n];
		} else {
			rcv_sample = rcv;
		}
	}
	bool* v_is_matched = new bool[2 * n];

	int* q_ptrs_start_sample = new int[n + 1];
	int* q_ptrs_end_sample = new int[n + 1];
	int* que_sample = new int[2 * n];
	int* que_locs_sample = new int[2 * n];
	int* degs_sample = new int[2 * n];

	int* unmatched_nbrs = new int[n];
	RT* unmatched_terms = new RT[n];

	const int adaptive_scaling_window = 100;
	const int adaptive_scaling_threshold = 99;
	const int adaptive_min_inner_scaling_iter = 0;
	const int adaptive_max_inner_scaling_iter = 10;
	int adaptive_current_inner_scaling_iter = no_inner_scaling_iter;
	int adaptive_window_samples = 0;
	int adaptive_window_perfect_matchings = 0;
	long long adaptive_increase_count = 0;
	long long adaptive_decrease_count = 0;

	PT estimation = 0;
	auto start_timex = high_resolution_clock::now();
	int x = 0; //sample ID
	while (true) {//get a sample
		x++;

		if constexpr(scale_approach > 1) {
			memcpy(rcv_sample, rcv, sizeof(ST) * 2 * n);
		}
		memset(v_is_matched, 0, sizeof(bool) * 2 * n);
		memcpy(q_ptrs_start_sample, q_ptrs_start, sizeof(int) * (n + 1));
		memcpy(q_ptrs_end_sample, q_ptrs_end, sizeof(int) * (n + 1));
		memcpy(que_sample, que, sizeof(int) * (2 * n));
		memcpy(que_locs_sample, que_locs, sizeof(int) * (2 * n));
		memcpy(degs_sample, degs, sizeof(int) * (2 * n));

		int min_deg_sample = min_deg;

		PT prod = 1; //estimation for sample - rv value

		int current_cardinality = 0;
		bool process_deferred_degree_2 = false;
		while(current_cardinality < n) {
			if(min_deg_sample > n) {
				if constexpr(defer_degree_2) {
					if(!process_deferred_degree_2) {
						process_deferred_degree_2 = true;
						min_deg_sample = 2;
						continue;
					}
				}
				break;
			}
			else if(q_ptrs_start_sample[min_deg_sample] == q_ptrs_end_sample[min_deg_sample]) {
				min_deg_sample++;
				continue;
			}

			if constexpr(defer_degree_2) {
				if(!process_deferred_degree_2 && min_deg_sample == 2) {
					min_deg_sample = 3;
					continue;
				}
			}

			int current_vertex = md_take_vertex_from_bucket_weighted(min_deg_sample,
																	q_ptrs_start_sample,
																	q_ptrs_end_sample,
																	que_sample,
																	que_locs_sample,
																	degs_sample,
																	v_is_matched,
																	ptrs,
																	ids,
																	A,
																	At,
																	nnz);
			if(current_vertex < 0) {
				min_deg_sample++;
				continue;
			}
			if(degs_sample[current_vertex] == 0) {
				continue;
			}

			//inner scale logic
			if constexpr(scale_approach == 2) {
				scale<ST, true, true, false, false>(A, At, rcv_sample, rcv_sample + n, v_is_matched, v_is_matched + n, 0, 0, no_inner_scaling_iter);
			} else if constexpr(scale_approach == 3) {
				light_scale<ST>(A, At, rcv_sample, rcv_sample + n, v_is_matched, v_is_matched + n, current_vertex);
			} else if constexpr(scale_approach == 4) {
				if(current_cardinality % hybrid_full_scale_period == 0) {
					scale<ST, true, true, false, false>(A, At, rcv_sample, rcv_sample + n, v_is_matched, v_is_matched + n, 0, 0, no_inner_scaling_iter);
				} else {
					light_scale<ST>(A, At, rcv_sample, rcv_sample + n, v_is_matched, v_is_matched + n, current_vertex);
				}
			} else if constexpr(scale_approach == 5) {
				scale<ST, true, true, false, false>(A, At, rcv_sample, rcv_sample + n, v_is_matched, v_is_matched + n, 0, 0, (2 * (double(current_cardinality)/n) * no_inner_scaling_iter) + 1);
			} else if constexpr(scale_approach == 6) {
				scale<ST, true, true, false, false>(A, At, rcv_sample, rcv_sample + n, v_is_matched, v_is_matched + n, 0, 0, (2  * no_inner_scaling_iter) -  (2 * (double(current_cardinality)/n) * no_inner_scaling_iter) + 1);
			} else if constexpr(scale_approach == 7) {
				if(adaptive_current_inner_scaling_iter > 0) {
					scale<ST, true, true, false, false>(A, At, rcv_sample, rcv_sample + n, v_is_matched, v_is_matched + n, 0, 0, adaptive_current_inner_scaling_iter);
				}
			} 

			int no_unmatched = 0; //find number of unmatched columns
			RT sum_unmatched = 0;
			for(int ptr = ptrs[current_vertex]; ptr < ptrs[current_vertex + 1]; ptr++) {
				int col_id = ids[ptr];
				if(v_is_matched[col_id] == false) {
					RT edge = (current_vertex < n)
						? static_cast<RT>(A.vals[ptr])
						: static_cast<RT>(At.vals[ptr - nnz]);
					RT proposal_term = edge;
					if constexpr(scale_approach > 0) {
						if constexpr(amplify) {
							proposal_term *= rcv_sample[col_id] * rcv_sample[col_id];
						} else {
							proposal_term *= rcv_sample[col_id];
						}
					}
					sum_unmatched += proposal_term;
					unmatched_terms[no_unmatched] = proposal_term;
					unmatched_nbrs[no_unmatched++] = col_id;
				}
			}
			if(no_unmatched == 0) {
				prod = 0;
				break; //if everything is matched sample contribution is 0
			}

			int chosen_column = -1;
			RT rand_val = sum_unmatched * dis(gen);
			chosen_column = unmatched_nbrs[0]; 
			RT running_sum = 0;
			for(int chosen = 0; chosen < no_unmatched; chosen++) {
				chosen_column = unmatched_nbrs[chosen];
				running_sum += unmatched_terms[chosen];
				if(rand_val <= running_sum) {
					break;
				}
			}
			if constexpr(scale_approach == 0) {
				prod = prod * sum_unmatched;
			} else if constexpr(amplify) {
				prod = prod * (sum_unmatched / (rcv_sample[chosen_column] * rcv_sample[chosen_column]));
			} else {
				prod = prod * (sum_unmatched / rcv_sample[chosen_column]);
			}

			v_is_matched[current_vertex] = v_is_matched[chosen_column] = true;
			degs_sample[current_vertex] = degs_sample[chosen_column] = 0;
			current_cardinality++;

			int delete_vertex = current_vertex;
			for(int rep = 0; rep < 2; rep++) {
				for(int ptr = ptrs[delete_vertex]; ptr < ptrs[delete_vertex + 1]; ptr++) {
					int nbr = ids[ptr];
					int deg = degs_sample[nbr];
					if(deg > 0) {
						int loc = que_locs_sample[nbr];
						int deg_head = que_sample[q_ptrs_start_sample[deg]++]; //extracts head shift the pointer
						
						que_sample[loc] = deg_head;  //put head to the loc
						que_locs_sample[deg_head] = loc; //set the locs of head
						
						que_sample[q_ptrs_end_sample[deg - 1]] = nbr; //put nbr to the end of the previous degree 
						que_locs_sample[nbr] = q_ptrs_end_sample[deg - 1]++; //shift the pointer

						degs_sample[nbr]--; //update the degree
						if(degs_sample[nbr] < min_deg_sample) {
							if constexpr(defer_degree_2) {
								if(process_deferred_degree_2 || degs_sample[nbr] != 2) {
									min_deg_sample = degs_sample[nbr];
								}
							} else {
								min_deg_sample = degs_sample[nbr];
							}
						}
					}
				}	
				delete_vertex = chosen_column;
			}
		}

		stats[current_cardinality]++;
		if(current_cardinality == n) {
			estimation += prod;
			products.push_back(prod);		
		} else {
			products.push_back(0);
		}

		if constexpr(scale_approach == 7) {
			adaptive_window_samples++;
			if(current_cardinality == n) {
				adaptive_window_perfect_matchings++;
			}

			if(adaptive_window_samples == adaptive_scaling_window) {
				if(adaptive_window_perfect_matchings < adaptive_scaling_threshold) {
					if(adaptive_current_inner_scaling_iter < adaptive_max_inner_scaling_iter) {
						adaptive_current_inner_scaling_iter++;
						adaptive_increase_count++;
					}
				} else {
					if(adaptive_current_inner_scaling_iter > adaptive_min_inner_scaling_iter) {
						adaptive_current_inner_scaling_iter--;
						adaptive_decrease_count++;
					}
				}
				adaptive_window_samples = 0;
				adaptive_window_perfect_matchings = 0;
			}
		}

		if(x % period == 0) {
			long double currentEstimation = estimation / x;
			auto duration = duration_cast<milliseconds>(high_resolution_clock::now() - start_timex);
			double seconds = ((double)duration.count()) / 1000;
			dataPoint p(x, currentEstimation, seconds);
			dataPoints.push_back(p);
			if(seconds > time_limit) {
				break;
			}
		}
	}
	auto stop_timex = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop_timex - start_timex);

	long int NO_SAMPLES = products.size();
	long double sum = 0;
	for(int i = 0; i < NO_SAMPLES; i++) {
		sum += products[i];
	}
	long double mean = sum / NO_SAMPLES;

	long double variance = 0;
	for(int i = 0; i < NO_SAMPLES; i++) {
		variance += (products[i] - mean) * (products[i] - mean);
	}
	variance /= (NO_SAMPLES - 1);
	long double stdev = sqrt(variance);

	long double stat_mean = 0;
	int q4pt = -1, q3pt = -1, q2pt = -1;
	long double sum_stat = 0;
	int max_cardinality = 0;
	for(long long i = 0; i <= n; i++) {
		sum_stat += stats[i];
		if((q4pt == -1) && ((sum_stat / NO_SAMPLES) > 0.25)) q4pt = i;
		if((q3pt == -1) && ((sum_stat / NO_SAMPLES) > 0.50)) q3pt = i;
		if((q2pt == -1) && ((sum_stat / NO_SAMPLES) > 0.75)) q2pt = i;
		if(stats[i] > 0) {
			max_cardinality = i;
		}
		stat_mean += i * stats[i];
	}
	stat_mean /= NO_SAMPLES;

	cout << "rasmussen_MD_"; 
	if constexpr(scale_approach == 0) {
		cout << "no_scaling";
	} else if constexpr(scale_approach == 1) {
		cout << "once_scaling";
	} else if constexpr(scale_approach == 2) {
		cout << "scaling_every_step";
	} else if constexpr(scale_approach == 3) {
		cout << "light_scaling_every_step";
	} else if constexpr(scale_approach == 4) {
		cout << "hybrid_scaling_every_step_period_" << hybrid_full_scale_period;
	} else if constexpr(scale_approach == 5) {
		cout << "scaling_every_step_with_interpolation";
	}  else if constexpr(scale_approach == 6) {
		cout << "scaling_every_step_with_rev_interpolation";
	}  else if constexpr(scale_approach == 7) {
		cout << "adaptive_scaling_every_step";
	}
	if constexpr(defer_degree_2) {
		cout << "_defer_degree_2";
	}
	if(md_bucket_scan_limit > 1) {
		cout << "_B" << md_bucket_scan_limit;
	}
	if(amplify) {
		cout << "_with_amplification";
	} 
	
	cout << " with " << NO_SAMPLES << " samples ------------------------------------------" << endl;
	cout << "Estimation: " << mean << "\tSamples_variance: " << variance << "\tSamples_stdev: " << stdev << endl;
	cout << "Ratio_perfect_matching: " << ((long double)stats[n]) / NO_SAMPLES << "\tMean_matching_cardinality: " << stat_mean << "\tQs: [" 
		 << q4pt << ", " << q3pt << ", " << q2pt << "] -- maxcard:" << max_cardinality << endl;
    cout << "Total_time_sec: " << ((double)duration.count()) / 1000 << "\tTime_taken_by_a_single_sample_ms: " << ((double)duration.count()) / NO_SAMPLES << endl << endl;
	if constexpr(scale_approach == 7) {
		cout << "Adaptive_scaling_window: " << adaptive_scaling_window
		     << "\tAdaptive_scaling_threshold: " << adaptive_scaling_threshold
		     << "\tAdaptive_final_inner_iterations: " << adaptive_current_inner_scaling_iter
		     << "\tAdaptive_increases: " << adaptive_increase_count
		     << "\tAdaptive_decreases: " << adaptive_decrease_count << endl << endl;
	}

	ofstream out_file;
	string d2_suffix = "";
	string run_bucket_suffix = md_bucket_suffix();
	if constexpr(defer_degree_2) {
		d2_suffix = "_deferD2";
	}
	if constexpr(scale_approach == 0) {
		if constexpr(amplify) {
			out_file.open (filename + ".MD_noscl" + d2_suffix + run_bucket_suffix + "_amplify.wout");
		} else {
	  		out_file.open (filename + ".MD_noscl" + d2_suffix + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 1) {
		if constexpr(amplify) {
			out_file.open (filename + ".MD_oncescl" + d2_suffix + run_bucket_suffix + "_amplify.wout");
		} else {
	  		out_file.open (filename + ".MD_oncescl" + d2_suffix + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 2) {
		if constexpr(amplify) {
			out_file.open (filename + ".MD_allscl" + d2_suffix + run_bucket_suffix + "_amplify.wout");
		} else {
  			out_file.open (filename + ".MD_allscl" + d2_suffix + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 3) {
		if constexpr(amplify) {
			out_file.open (filename + ".MD_lightscl" + d2_suffix + run_bucket_suffix + "_amplify.wout");
		} else {
  			out_file.open (filename + ".MD_lightscl" + d2_suffix + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 4) {
		if constexpr(amplify) {
			out_file.open (filename + ".MD_hybridscl_P" + to_string(hybrid_full_scale_period) + d2_suffix + run_bucket_suffix + "_amplify.wout");
		} else {
  			out_file.open (filename + ".MD_hybridscl_P" + to_string(hybrid_full_scale_period) + d2_suffix + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 5) {
		if constexpr(amplify) {
			out_file.open (filename + ".MD_intscl_P" + to_string(hybrid_full_scale_period) + d2_suffix + run_bucket_suffix + "_amplify.wout");
		} else {
  			out_file.open (filename + ".MD_intscl_P" + to_string(hybrid_full_scale_period) + d2_suffix + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 6) {
		if constexpr(amplify) {
			out_file.open (filename + ".MD_revintscl_P" + to_string(hybrid_full_scale_period) + d2_suffix + run_bucket_suffix + "_amplify.wout");
		} else {
  			out_file.open (filename + ".MD_revintscl_P" + to_string(hybrid_full_scale_period) + d2_suffix + run_bucket_suffix + ".wout");
		}
	} else if constexpr(scale_approach == 7) {
		if constexpr(amplify) {
			out_file.open (filename + ".MD_adaptivescl_W100_T99" + d2_suffix + run_bucket_suffix + "_amplify.wout");
		} else {
  			out_file.open (filename + ".MD_adaptivescl_W100_T99" + d2_suffix + run_bucket_suffix + ".wout");
		}
	}				
	out_file << n << "\n";
	for(int i = 1; i < n; i++) {
		out_file << stats[i] << " ";
	}
	out_file << stats[n] << "\n";
	out_file << dataPoints.size() << "\n";
	for(auto dp : dataPoints) {
		out_file << dp.noSamples << "\t" << dp.estimation << "\t" << dp.seconds << "\n";
	}
	out_file.close();

	delete [] ptrs;
	delete [] ids; 
	delete [] degs;
	delete [] degs_sample;
	delete [] deg_counts;
	delete [] que;
	delete [] que_sample;
	delete [] q_ptrs_start;
	delete [] q_ptrs_start_sample;
	delete [] q_ptrs_end;
	delete [] q_ptrs_end_sample;
	delete [] que_locs;
	delete [] que_locs_sample;
	delete [] stats;
	delete [] unmatched_nbrs;
	delete [] unmatched_terms;
	delete [] v_is_matched;

	if constexpr(scale_approach > 0) {
		delete [] rcv;
		if constexpr(scale_approach > 1) {
			delete [] rcv_sample;
		}
	}
}

int main(int argc, char *argv[]) {
	if(argc > 2) {
		filename = argv[1];
		time_limit = atof(argv[2]);
		if(argc > 3) {
			md_bucket_scan_limit = atoi(argv[3]);
		}
		if(argc > 4) {
			use_fixed_seed = true;
			fixed_seed = static_cast<unsigned int>(strtoul(argv[4], nullptr, 10));
		}
	} else {
		cout << "Usage: " << argv[0] << " <filename> <time_limit> [MD_bucket_scan_B] [seed]" << endl;
		exit(1);
	}
	if(md_bucket_scan_limit <= 0) {
		cout << "MD_bucket_scan_B must be > 0" << endl;
		exit(1);
	}

	Matrix<long double> matrix;
	read_coords(filename, matrix, n);
	print_matrix_stats(matrix, "Before tassa");

	matrix = tassa(matrix);
	if (matrix.empty()) {
		cout << "Permanent is 0 - no perfect matching" << endl;
		return 0;
	}
	print_matrix_stats(matrix, "After tassa");


	CSRMatrix A;
	convertDenseToSparse(matrix, A, n);
	CSRMatrix At = build_transpose(A);

	cout << "Matrix_dimension: " << A.n << " " << "\tMatrix_nnz: " << A.row_ptr[n] 
	     << "\tAvgNnz: " << ((double)A.row_ptr[n]) / n << endl << endl;
	cout << "MD_bucket_scan_B: " << md_bucket_scan_limit << endl << endl;
	if(use_fixed_seed) {
		cout << "RNG_seed: " << fixed_seed << endl << endl;
	} else {
		cout << "RNG_seed: random_device" << endl << endl;
	}


		/*//scale_approach is:
			//0: no scaling
			//1: once at the beginning
			//2: once + scaling at every step
			//3: once + light scaling at every step
			//4: once + hybrid scaling at every step
			//5: once + scaling at every step with interpolation
			//6: once + scaling at every step with reverse interpolation
			//7: once + adaptive scaling at every step
			template<typename PT, typename RT, typename ST, int scale_approach, bool amplify>
			void rasmussen_generic(CSRMatrix& A, CSRMatrix& At,
								int no_outer_scaling_iter, int no_inner_scaling_iter, 
								int hybrid_full_scale_period) {
		*/

	rasmussen_generic<long double, long double, long double, 0, false>(A, At, 0, 0, 0); //no scaling
	rasmussen_generic<long double, long double, long double, 1, false>(A, At, 50, 0, 0); //once scaling
	rasmussen_generic<long double, long double, long double, 2, false>(A, At, 50, 5, 0); //scaling every step
	rasmussen_generic<long double, long double, long double, 3, false>(A, At, 50, 5, 0); //light scaling at every step
	rasmussen_md_generic<long double, long double, long double, 0, false, false>(A, At, 0, 0, 0); //minimum degree without scaling
	rasmussen_md_generic<long double, long double, long double, 2, false, false>(A, At, 50, 5, 0); //scaling every step
	rasmussen_md_generic<long double, long double, long double, 3, false, false>(A, At, 50, 5, 0); //light scaling at every step
	rasmussen_md_generic<long double, long double, long double, 4, false, false>(A, At, 50, 5, 5); //hybrid scaling every step 5
	rasmussen_md_generic<long double, long double, long double, 4, false, false>(A, At, 50, 5, 10); //hybrid scaling every step 10 
	rasmussen_md_generic<long double, long double, long double, 4, false, false>(A, At, 50, 5, 20); //hybrid scaling every step 20

	return 0;
}
