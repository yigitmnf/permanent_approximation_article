#include <bits/stdc++.h>
using namespace std;
using Clock = chrono::steady_clock;

static double seconds_between(Clock::time_point a, Clock::time_point b) {
    return chrono::duration_cast<chrono::duration<double>>(b - a).count();
}

struct DenseMatrix {
    int n;
    vector<vector<long double>> a;
};

struct CsrMatrix {
    int n;
    vector<int> row_ptr;
    vector<int> col_idx;
    vector<long double> vals;
};

static DenseMatrix read_matrix_market(const string &path, bool force_pattern) {
    ifstream in(path);
    if (!in) { cerr << "dosya acilamadi: " << path << endl; exit(1); }
    bool symmetric_storage = false;
    bool pattern_storage = false;
    string line;
    while (getline(in, line)) {
        if (!line.empty() && line.rfind("%%MatrixMarket", 0) == 0) {
            if (line.find("symmetric") != string::npos) symmetric_storage = true;
            if (line.find("pattern") != string::npos) pattern_storage = true;
        }
        if (!line.empty() && line[0] != '%') break;
    }
    istringstream header(line);
    int rows, cols; long long entry_count;
    header >> rows >> cols >> entry_count;
    if (rows != cols) { cerr << "matris kare degil: " << path << endl; exit(1); }
    DenseMatrix m;
    m.n = rows;
    m.a.assign(rows, vector<long double>(rows, 0.0L));
    for (long long k = 0; k < entry_count; k++) {
        int i, j; long double value = 1.0L;
        if (pattern_storage) in >> i >> j;
        else { in >> i >> j >> value; }
        if (force_pattern) value = 1.0L;
        m.a[i - 1][j - 1] += value;
        if (symmetric_storage && i != j) m.a[j - 1][i - 1] += value;
    }
    if (force_pattern) {
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < rows; j++)
                if (m.a[i][j] != 0.0L) m.a[i][j] = 1.0L;
    }
    return m;
}

static vector<int> hopcroft_karp_dense(const DenseMatrix &m) {
    int n = m.n;
    const long double eps = 1e-9L;
    vector<vector<int>> adjacency(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (fabsl(m.a[i][j]) > eps) adjacency[i].push_back(j);
    vector<int> match_row(n, -1), match_col(n, -1), level(n);
    auto bfs_layers = [&]() {
        deque<int> queue_rows;
        bool reachable_free_col = false;
        for (int r = 0; r < n; r++) {
            if (match_row[r] == -1) { level[r] = 0; queue_rows.push_back(r); }
            else level[r] = INT_MAX;
        }
        while (!queue_rows.empty()) {
            int r = queue_rows.front(); queue_rows.pop_front();
            for (int c : adjacency[r]) {
                int r2 = match_col[c];
                if (r2 == -1) reachable_free_col = true;
                else if (level[r2] == INT_MAX) { level[r2] = level[r] + 1; queue_rows.push_back(r2); }
            }
        }
        return reachable_free_col;
    };
    function<bool(int)> augment = [&](int r) {
        for (int c : adjacency[r]) {
            int r2 = match_col[c];
            if (r2 == -1 || (level[r2] == level[r] + 1 && augment(r2))) {
                match_row[r] = c; match_col[c] = r; return true;
            }
        }
        level[r] = INT_MAX;
        return false;
    };
    while (bfs_layers())
        for (int r = 0; r < n; r++)
            if (match_row[r] == -1) augment(r);
    return match_row;
}

static vector<int> tarjan_scc_dense(const DenseMatrix &m) {
    int n = m.n;
    const long double eps = 1e-9L;
    vector<int> component(n, -1), index_of(n, -1), lowlink(n, 0), edge_iter(n, 0);
    vector<char> on_stack(n, 0);
    vector<int> stack_nodes, call_stack;
    int next_index = 0, component_count = 0;
    for (int start = 0; start < n; start++) {
        if (index_of[start] != -1) continue;
        call_stack.push_back(start);
        while (!call_stack.empty()) {
            int v = call_stack.back();
            if (index_of[v] == -1) {
                index_of[v] = lowlink[v] = next_index++;
                stack_nodes.push_back(v);
                on_stack[v] = 1;
                edge_iter[v] = 0;
            }
            bool descended = false;
            while (edge_iter[v] < n) {
                int w = edge_iter[v]++;
                if (w == v || fabsl(m.a[v][w]) <= eps) continue;
                if (index_of[w] == -1) { call_stack.push_back(w); descended = true; break; }
                if (on_stack[w]) lowlink[v] = min(lowlink[v], index_of[w]);
            }
            if (descended) continue;
            if (lowlink[v] == index_of[v]) {
                while (true) {
                    int w = stack_nodes.back(); stack_nodes.pop_back();
                    on_stack[w] = 0;
                    component[w] = component_count;
                    if (w == v) break;
                }
                component_count++;
            }
            call_stack.pop_back();
            if (!call_stack.empty()) {
                int parent = call_stack.back();
                lowlink[parent] = min(lowlink[parent], lowlink[v]);
            }
        }
    }
    return component;
}

static bool tassa_static_filter(DenseMatrix &m) {
    int n = m.n;
    vector<int> match_row = hopcroft_karp_dense(m);
    for (int i = 0; i < n; i++) if (match_row[i] == -1) return false;
    vector<int> row_of_col(n);
    for (int i = 0; i < n; i++) row_of_col[match_row[i]] = i;
    DenseMatrix permuted; permuted.n = n; permuted.a.resize(n);
    for (int i = 0; i < n; i++) permuted.a[i] = m.a[row_of_col[i]];
    vector<long double> diagonal_values(n);
    for (int i = 0; i < n; i++) { diagonal_values[i] = permuted.a[i][i]; permuted.a[i][i] = 0.0L; }
    vector<int> component = tarjan_scc_dense(permuted);
    DenseMatrix support; support.n = n;
    support.a.assign(n, vector<long double>(n, 0.0L));
    for (int i = 0; i < n; i++) support.a[i][i] = diagonal_values[i];
    const long double eps = 1e-9L;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (fabsl(permuted.a[i][j]) > eps && component[i] == component[j])
                support.a[i][j] = permuted.a[i][j];
    for (int i = 0; i < n; i++) m.a[match_row[i]] = support.a[i];
    return true;
}

static CsrMatrix dense_to_csr(const DenseMatrix &m) {
    const long double eps = 1e-9L;
    CsrMatrix c; c.n = m.n;
    c.row_ptr.assign(m.n + 1, 0);
    for (int i = 0; i < m.n; i++) {
        c.row_ptr[i] = (int)c.col_idx.size();
        for (int j = 0; j < m.n; j++)
            if (fabsl(m.a[i][j]) > eps) { c.col_idx.push_back(j); c.vals.push_back(m.a[i][j]); }
    }
    c.row_ptr[m.n] = (int)c.col_idx.size();
    return c;
}

static CsrMatrix transpose_csr(const CsrMatrix &a) {
    CsrMatrix t; t.n = a.n;
    t.row_ptr.assign(a.n + 1, 0);
    int nnz = a.row_ptr[a.n];
    t.col_idx.resize(nnz); t.vals.resize(nnz);
    for (int k = 0; k < nnz; k++) t.row_ptr[a.col_idx[k] + 1]++;
    for (int i = 0; i < a.n; i++) t.row_ptr[i + 1] += t.row_ptr[i];
    vector<int> cursor(t.row_ptr.begin(), t.row_ptr.end() - 1);
    for (int r = 0; r < a.n; r++)
        for (int p = a.row_ptr[r]; p < a.row_ptr[r + 1]; p++) {
            int c = a.col_idx[p];
            t.col_idx[cursor[c]] = r;
            t.vals[cursor[c]] = a.vals[p];
            cursor[c]++;
        }
    return t;
}

static void sinkhorn_scale(const CsrMatrix &A, const CsrMatrix &At,
                           long double *rv, long double *cv,
                           const bool *row_matched, const bool *col_matched,
                           int start_row, int iterations, bool initialize) {
    int n = A.n;
    if (initialize) for (int i = 0; i < n; i++) rv[i] = cv[i] = 1.0L;
    for (int it = 0; it < iterations; it++) {
        for (int c = 0; c < n; c++) {
            if (col_matched != nullptr && col_matched[c]) continue;
            long double sum = 0.0L;
            for (int p = At.row_ptr[c]; p < At.row_ptr[c + 1]; p++) {
                int r = At.col_idx[p];
                if (r < start_row) continue;
                if (row_matched != nullptr && row_matched[r]) continue;
                sum += At.vals[p] * rv[r];
            }
            if (sum != 0.0L) cv[c] = 1.0L / sum;
        }
        for (int r = start_row; r < n; r++) {
            if (row_matched != nullptr && row_matched[r]) continue;
            long double sum = 0.0L;
            for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; p++) {
                int c = A.col_idx[p];
                if (col_matched != nullptr && col_matched[c]) continue;
                sum += A.vals[p] * cv[c];
            }
            if (sum != 0.0L) rv[r] = 1.0L / sum;
        }
    }
}

static void light_scale(const CsrMatrix &A, const CsrMatrix &At,
                        long double *rv, long double *cv,
                        const bool *row_matched, const bool *col_matched,
                        int current_vertex) {
    int n = A.n;
    if (current_vertex < n) {
        int row = current_vertex;
        long double row_sum = 0.0L;
        for (int p = A.row_ptr[row]; p < A.row_ptr[row + 1]; p++) {
            int c = A.col_idx[p];
            if (col_matched[c]) continue;
            long double edge = A.vals[p];
            long double col_sum = 0.0L;
            for (int p2 = At.row_ptr[c]; p2 < At.row_ptr[c + 1]; p2++) {
                int r = At.col_idx[p2];
                if (row_matched[r]) continue;
                long double row_sum2 = 0.0L;
                for (int p3 = A.row_ptr[r]; p3 < A.row_ptr[r + 1]; p3++) {
                    int c2 = A.col_idx[p3];
                    if (col_matched[c2]) continue;
                    row_sum2 += A.vals[p3] * cv[c2];
                }
                if (row_sum2 != 0.0L) rv[r] = 1.0L / row_sum2;
                col_sum += At.vals[p2] * rv[r];
            }
            if (col_sum != 0.0L) cv[c] = 1.0L / col_sum;
            row_sum += edge * cv[c];
        }
        if (row_sum != 0.0L) rv[row] = 1.0L / row_sum;
    } else {
        int col = current_vertex - n;
        long double col_sum = 0.0L;
        for (int p = At.row_ptr[col]; p < At.row_ptr[col + 1]; p++) {
            int r = At.col_idx[p];
            if (row_matched[r]) continue;
            long double edge = At.vals[p];
            long double row_sum = 0.0L;
            for (int p2 = A.row_ptr[r]; p2 < A.row_ptr[r + 1]; p2++) {
                int c = A.col_idx[p2];
                if (col_matched[c]) continue;
                long double col_sum2 = 0.0L;
                for (int p3 = At.row_ptr[c]; p3 < At.row_ptr[c + 1]; p3++) {
                    int r2 = At.col_idx[p3];
                    if (row_matched[r2]) continue;
                    col_sum2 += At.vals[p3] * rv[r2];
                }
                if (col_sum2 != 0.0L) cv[c] = 1.0L / col_sum2;
                row_sum += A.vals[p2] * cv[c];
            }
            if (row_sum != 0.0L) rv[r] = 1.0L / row_sum;
            col_sum += edge * rv[r];
        }
        if (col_sum != 0.0L) cv[col] = 1.0L / col_sum;
    }
}

struct ShadowMatching {
    int n;
    const CsrMatrix *A;
    const CsrMatrix *At;
    vector<int> match_row, match_col;
    vector<int> visit_stamp;
    int current_stamp = 0;
    const vector<char> *edge_dead_row_side = nullptr;
    const bool *row_matched = nullptr;
    const bool *col_matched = nullptr;

    void init(const CsrMatrix &a, const CsrMatrix &at, const vector<int> &initial_match_row) {
        n = a.n; A = &a; At = &at;
        match_row = initial_match_row;
        match_col.assign(n, -1);
        for (int r = 0; r < n; r++) match_col[match_row[r]] = r;
        visit_stamp.assign(n, 0);
    }

    bool try_augment(int start_row) {
        current_stamp++;
        vector<pair<int,int>> dfs_stack;
        vector<int> parent_edge_row, parent_edge_col;
        dfs_stack.push_back({start_row, (*A).row_ptr[start_row]});
        vector<pair<int,int>> path;
        while (!dfs_stack.empty()) {
            auto &top = dfs_stack.back();
            int r = top.first;
            bool advanced = false;
            for (int &p = top.second; p < (*A).row_ptr[r + 1]; ) {
                int c = (*A).col_idx[p];
                int edge_index = p;
                p++;
                if (col_matched[c]) continue;
                if (edge_dead_row_side != nullptr && (*edge_dead_row_side)[edge_index]) continue;
                if (visit_stamp[c] == current_stamp) continue;
                visit_stamp[c] = current_stamp;
                int r2 = match_col[c];
                path.push_back({r, c});
                if (r2 == -1) {
                    for (auto &pr : path) { match_row[pr.first] = pr.second; match_col[pr.second] = pr.first; }
                    return true;
                }
                dfs_stack.push_back({r2, (*A).row_ptr[r2]});
                advanced = true;
                break;
            }
            if (!advanced) {
                dfs_stack.pop_back();
                if (!path.empty()) path.pop_back();
            }
        }
        return false;
    }

    bool remove_pair_and_repair(int sampled_row, int sampled_col) {
        int paired_col = match_row[sampled_row];
        int paired_row = match_col[sampled_col];
        if (paired_col == sampled_col) {
            match_row[sampled_row] = -1;
            match_col[sampled_col] = -1;
            return true;
        }
        match_row[sampled_row] = -1; match_col[paired_col] = -1;
        match_row[paired_row] = -1; match_col[sampled_col] = -1;
        return try_augment(paired_row);
    }
};

struct RunOptions {
    string matrix_path;
    string variant;
    long long target_samples = 10000;
    double time_cap_seconds = 600.0;
    unsigned seed = 42;
    string output_prefix;
    bool force_pattern = false;
    bool record_curve = false;
    bool record_doom = false;
    string tassa_mode = "never";
    long double tassa_threshold = 1e30L;
    bool record_ds = false;
    int outer_scaling_iterations = 50;
    int inner_scaling_iterations = 5;
};

struct RunResult {
    long long samples = 0, perfect_matchings = 0;
    long double sum_cardinality = 0;
    vector<long long> cardinality_histogram;
    long double best_log_product = -INFINITY;
    long double exp_shift_sum = 0;
    double total_seconds = 0, scaling_seconds = 0, tassa_seconds = 0;
    long long tassa_calls = 0, tassa_kills = 0, tassa_early_zero = 0, tassa_trigger_opportunities = 0;
    vector<double> curve_alive_sum;
    vector<double> ds_dev_sum;
    vector<long long> ds_dev_cnt;
    vector<array<double,3>> doom_records;
    long long doom_born_count = 0;

    void add_log_product(long double log10_product) {
        if (log10_product > best_log_product) {
            exp_shift_sum = exp_shift_sum * powl(10.0L, best_log_product - log10_product) + 1.0L;
            best_log_product = log10_product;
        } else {
            exp_shift_sum += powl(10.0L, log10_product - best_log_product);
        }
    }

    long double estimate_log10() const {
        if (perfect_matchings == 0 || samples == 0) return -INFINITY;
        return best_log_product + log10l(exp_shift_sum) - log10l((long double)samples);
    }
};

static void run_row_order_sampler(const CsrMatrix &A, const CsrMatrix &At, int scale_mode,
                                  const RunOptions &opt, const vector<int> &initial_matching,
                                  RunResult &out) {
    int n = A.n;
    mt19937 rng(opt.seed);
    uniform_real_distribution<double> uniform01(0.0, 1.0);
    vector<char> col_matched_flags(n), row_matched_flags(n);
    vector<int> candidate_columns(n);
    vector<long double> candidate_terms(n);
    vector<long double> base_rv(n, 1.0L), base_cv(n, 1.0L);
    vector<long double> rv(n), cv(n);
    int nnz = A.row_ptr[n];

    auto scaling_start = Clock::now();
    if (scale_mode > 0)
        sinkhorn_scale(A, At, base_rv.data(), base_cv.data(), nullptr, nullptr, 0, opt.outer_scaling_iterations, true);
    out.scaling_seconds += seconds_between(scaling_start, Clock::now());

    ShadowMatching shadow;
    bool shadow_active = opt.record_doom;
    out.cardinality_histogram.assign(n + 1, 0);
    if (opt.record_curve) out.curve_alive_sum.assign(n + 1, 0.0);

    auto run_start = Clock::now();
    for (long long sample_id = 1; sample_id <= opt.target_samples; sample_id++) {
        fill(col_matched_flags.begin(), col_matched_flags.end(), 0);
        fill(row_matched_flags.begin(), row_matched_flags.end(), 0);
        if (scale_mode >= 2) { rv = base_rv; cv = base_cv; }
        long double *rv_use = (scale_mode >= 2) ? rv.data() : base_rv.data();
        long double *cv_use = (scale_mode >= 2) ? cv.data() : base_cv.data();
        if (shadow_active) {
            shadow.init(A, At, initial_matching);
            shadow.row_matched = (bool*)row_matched_flags.data();
            shadow.col_matched = (bool*)col_matched_flags.data();
        }
        long long alive_edges = nnz;
        if (opt.record_curve) out.curve_alive_sum[0] += (double)alive_edges;

        long double log10_product = 0.0L;
        int doom_birth_iteration = -1;
        long double doom_edge_scaling_value = -1.0L;
        int completed = 0;
        int step;
        for (step = 0; step < n; step++) {
            if (scale_mode == 2) {
                auto t0 = Clock::now();
                sinkhorn_scale(A, At, rv_use, cv_use, nullptr, (bool*)col_matched_flags.data(), step, opt.inner_scaling_iterations, false);
                out.scaling_seconds += seconds_between(t0, Clock::now());
            } else if (scale_mode == 3) {
                auto t0 = Clock::now();
                light_scale(A, At, rv_use, cv_use, (bool*)row_matched_flags.data(), (bool*)col_matched_flags.data(), step);
                out.scaling_seconds += seconds_between(t0, Clock::now());
            }
            int candidate_count = 0;
            long double term_sum = 0.0L, term_max = 0.0L;
            for (int p = A.row_ptr[step]; p < A.row_ptr[step + 1]; p++) {
                int c = A.col_idx[p];
                if (col_matched_flags[c]) continue;
                long double term = A.vals[p];
                if (scale_mode > 0) term *= cv_use[c];
                candidate_columns[candidate_count] = c;
                candidate_terms[candidate_count++] = term;
                term_sum += term;
                if (term > term_max) term_max = term;
            }
            if (candidate_count == 0) break;
            long double draw = term_sum * (long double)uniform01(rng);
            int chosen_column = candidate_columns[0];
            long double running = 0.0L, chosen_term = candidate_terms[0];
            for (int k = 0; k < candidate_count; k++) {
                chosen_column = candidate_columns[k];
                chosen_term = candidate_terms[k];
                running += candidate_terms[k];
                if (draw <= running) break;
            }
            if (scale_mode == 0) log10_product += log10l(term_sum);
            else log10_product += log10l(term_sum / cv_use[chosen_column]);

            row_matched_flags[step] = 1;
            col_matched_flags[chosen_column] = 1;
            alive_edges -= candidate_count;
            for (int p = At.row_ptr[chosen_column]; p < At.row_ptr[chosen_column + 1]; p++) {
                int r = At.col_idx[p];
                if (r > step) alive_edges--;
            }
            if (opt.record_curve) out.curve_alive_sum[step + 1] += (double)alive_edges;

            if (shadow_active && doom_birth_iteration < 0) {
                if (!shadow.remove_pair_and_repair(step, chosen_column)) {
                    doom_birth_iteration = step;
                    long double raw_a = 0.0L;
                    for (int p = A.row_ptr[step]; p < A.row_ptr[step + 1]; p++)
                        if (A.col_idx[p] == chosen_column) raw_a = A.vals[p];
                    doom_edge_scaling_value = (scale_mode > 0)
                        ? rv_use[step] * raw_a * cv_use[chosen_column]
                        : -1.0L;
                }
            }
            if (step == n - 1) completed = 1;
        }
        int cardinality = completed ? n : step;
        out.cardinality_histogram[cardinality]++;
        out.sum_cardinality += cardinality;
        out.samples++;
        if (completed) { out.perfect_matchings++; out.add_log_product(log10_product); }
        else if (opt.record_doom && doom_birth_iteration >= 0) {
            out.doom_born_count++;
            out.doom_records.push_back({(double)doom_birth_iteration, (double)cardinality, (double)doom_edge_scaling_value});
        }
        if (sample_id % 64 == 0 && seconds_between(run_start, Clock::now()) > opt.time_cap_seconds) break;
    }
    out.total_seconds = seconds_between(run_start, Clock::now());
}

struct MinDegreeState {
    int n, nnz;
    vector<int> ptrs, ids, edge_of_ptr;
    vector<int> base_degrees, base_queue, base_queue_locs, base_bucket_start, base_bucket_end;
    int base_min_degree;

    void build(const CsrMatrix &A, const CsrMatrix &At) {
        n = A.n; nnz = A.row_ptr[n];
        ptrs.assign(2 * n + 1, 0);
        ids.assign(2 * nnz, 0);
        edge_of_ptr.assign(2 * nnz, 0);
        for (int i = 0; i < n; i++) { ptrs[i] = A.row_ptr[i]; ptrs[i + n] = At.row_ptr[i] + nnz; }
        ptrs[2 * n] = 2 * nnz;
        for (int k = 0; k < nnz; k++) { ids[k] = A.col_idx[k] + n; edge_of_ptr[k] = k; }
        vector<int> cursor(n);
        for (int c = 0; c < n; c++) cursor[c] = At.row_ptr[c];
        for (int r = 0; r < n; r++)
            for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; p++) {
                int c = A.col_idx[p];
                int col_side_position = cursor[c]++;
                ids[nnz + col_side_position] = r;
                edge_of_ptr[nnz + col_side_position] = p;
            }
        base_degrees.assign(2 * n, 0);
        vector<int> degree_counts(n + 1, 0);
        for (int v = 0; v < 2 * n; v++) { base_degrees[v] = ptrs[v + 1] - ptrs[v]; degree_counts[base_degrees[v]]++; }
        base_bucket_start.assign(n + 1, 0);
        base_bucket_end.assign(n + 1, 0);
        for (int d = 1; d <= n; d++) base_bucket_start[d] = base_bucket_start[d - 1] + degree_counts[d - 1];
        for (int d = 0; d <= n; d++) base_bucket_end[d] = base_bucket_start[d] + degree_counts[d];
        base_queue.assign(2 * n, 0);
        base_queue_locs.assign(2 * n, 0);
        vector<int> fill_cursor(base_bucket_start);
        for (int v = 0; v < 2 * n; v++) {
            int loc = fill_cursor[base_degrees[v]]++;
            base_queue_locs[v] = loc;
            base_queue[loc] = v;
        }
        base_min_degree = 0;
        for (int d = 0; d <= n; d++) if (degree_counts[d] > 0) { base_min_degree = d; break; }
    }
};

static void run_min_degree_sampler(const CsrMatrix &A, const CsrMatrix &At, int scale_mode,
                                   const RunOptions &opt, const vector<int> &initial_matching,
                                   RunResult &out) {
    int n = A.n;
    mt19937 rng(opt.seed);
    uniform_real_distribution<double> uniform01(0.0, 1.0);
    MinDegreeState md;
    md.build(A, At);
    int nnz = md.nnz;

    vector<long double> base_rcv(2 * n, 1.0L), rcv(2 * n);
    auto scaling_start = Clock::now();
    if (scale_mode > 0)
        sinkhorn_scale(A, At, base_rcv.data(), base_rcv.data() + n, nullptr, nullptr, 0, opt.outer_scaling_iterations, true);
    out.scaling_seconds += seconds_between(scaling_start, Clock::now());

    vector<char> vertex_matched(2 * n);
    vector<int> queue_s(2 * n), queue_locs_s(2 * n), degrees_s(2 * n), bucket_start_s(n + 1), bucket_end_s(n + 1);
    vector<int> candidate_ids(n), candidate_edges(n);
    vector<long double> candidate_terms(n);
    vector<char> edge_dead(nnz);
    bool tassa_on = (opt.tassa_mode != "never");

    ShadowMatching shadow;
    bool shadow_active = opt.record_doom || tassa_on;

    vector<int> residual_component(n, -1);
    vector<int> scc_index(n), scc_lowlink(n), scc_stack, scc_call_stack, scc_edge_iter(n);
    vector<char> scc_on_stack(n);

    out.cardinality_histogram.assign(n + 1, 0);
    if (opt.record_curve) out.curve_alive_sum.assign(n + 1, 0.0);
    int ds_stride = max(1, n / 25);
    if (opt.record_ds) { out.ds_dev_sum.assign(n + 1, 0.0); out.ds_dev_cnt.assign(n + 1, 0); }

    auto bucket_decrement = [&](int vertex, int &min_degree_sample) {
        int d = degrees_s[vertex];
        int loc = queue_locs_s[vertex];
        int head = queue_s[bucket_start_s[d]++];
        queue_s[loc] = head;
        queue_locs_s[head] = loc;
        queue_s[bucket_end_s[d - 1]] = vertex;
        queue_locs_s[vertex] = bucket_end_s[d - 1]++;
        degrees_s[vertex]--;
        if (degrees_s[vertex] < min_degree_sample) min_degree_sample = degrees_s[vertex];
    };

    auto tassa_filter_residual = [&](int &min_degree_sample, long long &alive_edges) -> long long {
        auto t0 = Clock::now();
        for (int c = 0; c < n; c++)
            residual_component[c] = vertex_matched[n + c] ? -1 : -2;
        int next_index = 0;
        int component_count = 0;
        for (int c = 0; c < n; c++) scc_index[c] = -1;
        for (int start_col = 0; start_col < n; start_col++) {
            if (vertex_matched[n + start_col] || scc_index[start_col] != -1) continue;
            scc_call_stack.push_back(start_col);
            while (!scc_call_stack.empty()) {
                int c = scc_call_stack.back();
                if (scc_index[c] == -1) {
                    scc_index[c] = scc_lowlink[c] = next_index++;
                    scc_stack.push_back(c); scc_on_stack[c] = 1;
                    scc_edge_iter[c] = A.row_ptr[shadow.match_col[c]];
                }
                int owner_row = shadow.match_col[c];
                bool descended = false;
                while (scc_edge_iter[c] < A.row_ptr[owner_row + 1]) {
                    int p = scc_edge_iter[c]++;
                    int c2 = A.col_idx[p];
                    if (c2 == c || vertex_matched[n + c2] || edge_dead[p]) continue;
                    if (scc_index[c2] == -1) { scc_call_stack.push_back(c2); descended = true; break; }
                    if (scc_on_stack[c2]) scc_lowlink[c] = min(scc_lowlink[c], scc_index[c2]);
                }
                if (descended) continue;
                if (scc_lowlink[c] == scc_index[c]) {
                    while (true) {
                        int w = scc_stack.back(); scc_stack.pop_back();
                        scc_on_stack[w] = 0;
                        residual_component[w] = component_count;
                        if (w == c) break;
                    }
                    component_count++;
                }
                scc_call_stack.pop_back();
                if (!scc_call_stack.empty()) {
                    int parent = scc_call_stack.back();
                    scc_lowlink[parent] = min(scc_lowlink[parent], scc_lowlink[c]);
                }
            }
        }
        long long kills = 0;
        for (int r = 0; r < n; r++) {
            if (vertex_matched[r]) continue;
            int matched_column_of_row = shadow.match_row[r];
            for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; p++) {
                int c = A.col_idx[p];
                if (vertex_matched[n + c] || edge_dead[p]) continue;
                if (c == matched_column_of_row) continue;
                if (residual_component[matched_column_of_row] != residual_component[c]) {
                    edge_dead[p] = 1;
                    kills++;
                    alive_edges--;
                    bucket_decrement(r, min_degree_sample);
                    bucket_decrement(c + n, min_degree_sample);
                }
            }
        }
        out.tassa_seconds += seconds_between(t0, Clock::now());
        out.tassa_calls++;
        out.tassa_kills += kills;
        return kills;
    };

    auto run_start = Clock::now();
    for (long long sample_id = 1; sample_id <= opt.target_samples; sample_id++) {
        if (scale_mode >= 2) rcv = base_rcv;
        long double *rcv_use = (scale_mode >= 2) ? rcv.data() : base_rcv.data();
        fill(vertex_matched.begin(), vertex_matched.end(), 0);
        copy(md.base_bucket_start.begin(), md.base_bucket_start.end(), bucket_start_s.begin());
        copy(md.base_bucket_end.begin(), md.base_bucket_end.end(), bucket_end_s.begin());
        copy(md.base_queue.begin(), md.base_queue.end(), queue_s.begin());
        copy(md.base_queue_locs.begin(), md.base_queue_locs.end(), queue_locs_s.begin());
        copy(md.base_degrees.begin(), md.base_degrees.end(), degrees_s.begin());
        int min_degree_sample = md.base_min_degree;
        if (tassa_on) fill(edge_dead.begin(), edge_dead.end(), 0);
        if (shadow_active) {
            shadow.init(A, At, initial_matching);
            shadow.row_matched = (bool*)vertex_matched.data();
            shadow.col_matched = (bool*)(vertex_matched.data() + n);
            shadow.edge_dead_row_side = tassa_on ? &edge_dead : nullptr;
        }
        long long alive_edges = nnz;
        if (opt.record_curve) out.curve_alive_sum[0] += (double)alive_edges;

        long double log10_product = 0.0L;
        int cardinality = 0;
        int doom_birth_iteration = -1;
        long double doom_edge_scaling_value = -1.0L;
        bool early_zero = false;

        while (cardinality < n) {
            if (min_degree_sample > n) break;
            if (bucket_start_s[min_degree_sample] == bucket_end_s[min_degree_sample]) { min_degree_sample++; continue; }
            while (bucket_start_s[min_degree_sample] < bucket_end_s[min_degree_sample] &&
                   degrees_s[queue_s[bucket_start_s[min_degree_sample]]] == 0)
                bucket_start_s[min_degree_sample]++;
            if (bucket_start_s[min_degree_sample] == bucket_end_s[min_degree_sample]) { min_degree_sample++; continue; }
            int current_vertex = queue_s[bucket_start_s[min_degree_sample]++];
            if (degrees_s[current_vertex] == 0) continue;

            if (scale_mode == 2) {
                auto t0 = Clock::now();
                sinkhorn_scale(A, At, rcv_use, rcv_use + n, (bool*)vertex_matched.data(), (bool*)(vertex_matched.data() + n), 0, opt.inner_scaling_iterations, false);
                out.scaling_seconds += seconds_between(t0, Clock::now());
            } else if (scale_mode == 3) {
                auto t0 = Clock::now();
                light_scale(A, At, rcv_use, rcv_use + n, (bool*)vertex_matched.data(), (bool*)(vertex_matched.data() + n), current_vertex);
                out.scaling_seconds += seconds_between(t0, Clock::now());
            }

            int candidate_count = 0;
            long double term_sum = 0.0L, term_max = 0.0L;
            for (int p = md.ptrs[current_vertex]; p < md.ptrs[current_vertex + 1]; p++) {
                int neighbor = md.ids[p];
                if (vertex_matched[neighbor]) continue;
                int edge_index = md.edge_of_ptr[p];
                if (tassa_on && edge_dead[edge_index]) continue;
                long double edge_value = (current_vertex < n) ? A.vals[p] : At.vals[p - nnz];
                long double term = edge_value;
                if (scale_mode > 0) term *= rcv_use[neighbor];
                candidate_ids[candidate_count] = neighbor;
                candidate_edges[candidate_count] = edge_index;
                candidate_terms[candidate_count++] = term;
                term_sum += term;
                if (term > term_max) term_max = term;
            }
            if (candidate_count == 0) break;
            long double draw = term_sum * (long double)uniform01(rng);
            int chosen_neighbor = candidate_ids[0];
            long double running = 0.0L, chosen_term = candidate_terms[0];
            for (int k = 0; k < candidate_count; k++) {
                chosen_neighbor = candidate_ids[k];
                chosen_term = candidate_terms[k];
                running += candidate_terms[k];
                if (draw <= running) break;
            }
            if (scale_mode == 0) log10_product += log10l(term_sum);
            else log10_product += log10l(term_sum / rcv_use[chosen_neighbor]);

            int deg_current = degrees_s[current_vertex];
            int deg_chosen = degrees_s[chosen_neighbor];
            vertex_matched[current_vertex] = vertex_matched[chosen_neighbor] = 1;
            degrees_s[current_vertex] = degrees_s[chosen_neighbor] = 0;
            cardinality++;
            alive_edges -= (deg_current + deg_chosen - 1);

            int delete_vertex = current_vertex;
            for (int rep = 0; rep < 2; rep++) {
                for (int p = md.ptrs[delete_vertex]; p < md.ptrs[delete_vertex + 1]; p++) {
                    int neighbor = md.ids[p];
                    if (degrees_s[neighbor] > 0) {
                        if (tassa_on && edge_dead[md.edge_of_ptr[p]]) continue;
                        bucket_decrement(neighbor, min_degree_sample);
                    }
                }
                delete_vertex = chosen_neighbor;
            }
            if (opt.record_curve && cardinality <= n) out.curve_alive_sum[cardinality] += (double)alive_edges;
            if (opt.record_ds && scale_mode > 0 && cardinality <= n && cardinality % ds_stride == 0) {
                long double dev_total = 0.0L;
                int free_rows = 0;
                for (int i = 0; i < n; i++) {
                    if (vertex_matched[i]) continue;
                    long double rs = 0.0L;
                    for (int p = A.row_ptr[i]; p < A.row_ptr[i + 1]; p++) {
                        int j = A.col_idx[p];
                        if (vertex_matched[n + j]) continue;
                        rs += rcv_use[i] * A.vals[p] * rcv_use[n + j];
                    }
                    dev_total += fabsl(rs - 1.0L);
                    free_rows++;
                }
                if (free_rows > 0) {
                    out.ds_dev_sum[cardinality] += (double)(dev_total / (long double)free_rows);
                    out.ds_dev_cnt[cardinality]++;
                }
            }

            int sampled_row, sampled_col;
            if (current_vertex < n) { sampled_row = current_vertex; sampled_col = chosen_neighbor - n; }
            else { sampled_row = chosen_neighbor; sampled_col = current_vertex - n; }

            if (shadow_active) {
                if (doom_birth_iteration < 0) {
                    if (!shadow.remove_pair_and_repair(sampled_row, sampled_col)) {
                        doom_birth_iteration = cardinality - 1;
                        long double raw_a = 0.0L;
                        for (int p = A.row_ptr[sampled_row]; p < A.row_ptr[sampled_row + 1]; p++)
                            if (A.col_idx[p] == sampled_col) raw_a = A.vals[p];
                        doom_edge_scaling_value = (scale_mode > 0)
                            ? rcv_use[sampled_row] * raw_a * rcv_use[sampled_col + n]
                            : -1.0L;
                        if (tassa_on) { early_zero = true; out.tassa_early_zero++; }
                    }
                }
                if (early_zero) break;
                if (tassa_on && doom_birth_iteration < 0 && cardinality < n) {
                    bool trigger = false;
                    if (opt.tassa_mode == "always") trigger = true;
                    else if (opt.tassa_mode == "adaptive") {
                        out.tassa_trigger_opportunities++;
                        if (candidate_count >= 2 && chosen_term < term_max &&
                            rcv_use[current_vertex] * chosen_term < opt.tassa_threshold) trigger = true;
                    }
                    if (trigger) tassa_filter_residual(min_degree_sample, alive_edges);
                }
            }
        }
        out.cardinality_histogram[cardinality]++;
        out.sum_cardinality += cardinality;
        out.samples++;
        if (cardinality == n) out.add_log_product(log10_product);
        else if (opt.record_doom && doom_birth_iteration >= 0) {
            out.doom_born_count++;
            out.doom_records.push_back({(double)doom_birth_iteration, (double)cardinality, (double)doom_edge_scaling_value});
        }
        if (cardinality == n) out.perfect_matchings++;
        if (sample_id % 64 == 0 && seconds_between(run_start, Clock::now()) > opt.time_cap_seconds) break;
    }
    out.total_seconds = seconds_between(run_start, Clock::now());
}

int main(int argc, char **argv) {
    if (argc < 7) {
        cerr << "kullanim: " << argv[0]
             << " <mtx> <varyant> <ornek_sayisi> <sure_limiti_sn> <tohum> <cikti_oneki>"
             << " [--pattern] [--egri] [--dogum] [--tassa never|adaptive|always] [--esik tau] [--dsizle]" << endl;
        cerr << "varyant: ras once light every md md_once md_light md_every" << endl;
        return 1;
    }
    RunOptions opt;
    opt.matrix_path = argv[1];
    opt.variant = argv[2];
    opt.target_samples = atoll(argv[3]);
    opt.time_cap_seconds = atof(argv[4]);
    opt.seed = (unsigned)strtoul(argv[5], nullptr, 10);
    opt.output_prefix = argv[6];
    for (int i = 7; i < argc; i++) {
        string flag = argv[i];
        if (flag == "--pattern") opt.force_pattern = true;
        else if (flag == "--egri") opt.record_curve = true;
        else if (flag == "--dogum") opt.record_doom = true;
        else if (flag == "--tassa" && i + 1 < argc) opt.tassa_mode = argv[++i];
        else if (flag == "--esik" && i + 1 < argc) opt.tassa_threshold = strtold(argv[++i], nullptr);
        else if (flag == "--dsizle") opt.record_ds = true;
    }

    DenseMatrix dense = read_matrix_market(opt.matrix_path, opt.force_pattern);
    int nnz_before = 0;
    for (auto &row : dense.a) for (auto v : row) if (fabsl(v) > 1e-9L) nnz_before++;
    if (!tassa_static_filter(dense)) {
        ofstream summary(opt.output_prefix + ".ozet");
        summary << "matris\t" << opt.matrix_path << "\tvaryant\t" << opt.variant
                << "\tsonuc\tPERMANENT_SIFIR" << endl;
        return 0;
    }
    CsrMatrix A = dense_to_csr(dense);
    CsrMatrix At = transpose_csr(A);
    vector<int> initial_matching = hopcroft_karp_dense(dense);

    int scale_mode = 0;
    bool use_min_degree = false;
    if (opt.variant == "ras") scale_mode = 0;
    else if (opt.variant == "once") scale_mode = 1;
    else if (opt.variant == "every") scale_mode = 2;
    else if (opt.variant == "light") scale_mode = 3;
    else if (opt.variant == "md") { scale_mode = 0; use_min_degree = true; }
    else if (opt.variant == "md_once") { scale_mode = 1; use_min_degree = true; }
    else if (opt.variant == "md_every") { scale_mode = 2; use_min_degree = true; }
    else if (opt.variant == "md_light") { scale_mode = 3; use_min_degree = true; }
    else { cerr << "bilinmeyen varyant: " << opt.variant << endl; return 1; }

    RunResult result;
    if (use_min_degree) run_min_degree_sampler(A, At, scale_mode, opt, initial_matching, result);
    else run_row_order_sampler(A, At, scale_mode, opt, initial_matching, result);

    int n = A.n;
    long long cumulative = 0;
    int q25 = -1, q50 = -1, q75 = -1;
    for (int c = 0; c <= n; c++) {
        cumulative += result.cardinality_histogram[c];
        if (q25 == -1 && cumulative * 4 >= result.samples) q25 = c;
        if (q50 == -1 && cumulative * 2 >= result.samples) q50 = c;
        if (q75 == -1 && cumulative * 4 >= 3 * result.samples) q75 = c;
    }

    ofstream summary(opt.output_prefix + ".ozet");
    summary.setf(ios::fixed);
    summary << setprecision(9);
    summary << "matris\t" << opt.matrix_path << "\n";
    summary << "varyant\t" << opt.variant << "\n";
    summary << "tassa_modu\t" << opt.tassa_mode << "\n";
    summary << "n\t" << n << "\n";
    summary << "nnz_ham\t" << nnz_before << "\n";
    summary << "nnz_tassa_sonrasi\t" << A.row_ptr[n] << "\n";
    summary << "tohum\t" << opt.seed << "\n";
    summary << "ornek\t" << result.samples << "\n";
    summary << "tam_eslesme\t" << result.perfect_matchings << "\n";
    summary << "tam_eslesme_orani\t" << (double)result.perfect_matchings / max(1LL, result.samples) << "\n";
    summary << "ortalama_kardinalite\t" << (double)(result.sum_cardinality / max(1LL, result.samples)) << "\n";
    summary << "kardinalite_q25_q50_q75\t" << q25 << " " << q50 << " " << q75 << "\n";
    long double est = result.estimate_log10();
    if (isinf(est)) summary << "tahmin_log10\tYOK\n";
    else summary << "tahmin_log10\t" << (double)est << "\n";
    summary << "sure_toplam_sn\t" << result.total_seconds << "\n";
    summary << "sure_scaling_sn\t" << result.scaling_seconds << "\n";
    summary << "sure_tassa_sn\t" << result.tassa_seconds << "\n";
    summary << "tassa_cagri\t" << result.tassa_calls << "\n";
    summary << "tassa_silinen\t" << result.tassa_kills << "\n";
    summary << "tassa_erken_sifir\t" << result.tassa_early_zero << "\n";
    summary << "tassa_firsat\t" << result.tassa_trigger_opportunities << "\n";
    summary << "dogum_kaydi\t" << result.doom_born_count << "\n";
    summary << "ornek_basina_ms\t" << 1000.0 * result.total_seconds / max(1LL, result.samples) << "\n";
    summary.close();

    if (opt.record_curve) {
        ofstream curve(opt.output_prefix + ".egri");
        for (int t = 0; t <= n; t++)
            curve << t << "\t" << result.curve_alive_sum[t] / max(1LL, result.samples) << "\n";
    }
    if (opt.record_doom) {
        ofstream doom(opt.output_prefix + ".dogum");
        doom << setprecision(9);
        for (auto &rec : result.doom_records)
            doom << (long long)rec[0] << "\t" << (long long)rec[1] << "\t" << n << "\t" << (double)rec[2] << "\n";
    }
    if (opt.record_ds) {
        ofstream ds(opt.output_prefix + ".dsat");
        ds << setprecision(9);
        for (int t = 0; t < (int)result.ds_dev_cnt.size(); t++)
            if (result.ds_dev_cnt[t] > 0)
                ds << t << "\t" << result.ds_dev_sum[t] / (double)result.ds_dev_cnt[t] << "\t" << result.ds_dev_cnt[t] << "\n";
    }
    return 0;
}
