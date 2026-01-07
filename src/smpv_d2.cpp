#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

struct Triplet { int row, col; double val; };

struct CsrMatrix {
    int rows = 0;
    int cols = 0;
    int64_t nnz = 0;
    vector<int> row_ptr;
    vector<int> col_ind;
    vector<double> values;
};

static vector<int> g_col_map;

static inline void require(bool cond, const string& msg) {
    if (!cond) throw runtime_error(msg);
}

static inline int owner_of(int global_idx, int P) {
    int o = global_idx % P;
    return (o < 0) ? (o + P) : o;
}

static inline int local_index_mod(int global_idx, int rank, int P) {
    return (global_idx - rank) / P;
}

static inline int local_count_mod(int global_len, int rank, int P) {
    if (global_len <= 0) return 0;
    int base = global_len / P;
    int rem  = global_len % P;
    return base + ((rank < rem) ? 1 : 0);
}

static inline string basename_noext(const string& path) {
    size_t pos = path.find_last_of("/\\");
    string fn = (pos == string::npos) ? path : path.substr(pos + 1);
    const string ext = ".mtx";
    if (fn.size() > ext.size() && fn.compare(fn.size() - ext.size(), ext.size(), ext) == 0)
        fn.erase(fn.size() - ext.size(), ext.size());
    return fn;
}

static inline int omp_threads_active() {
    int t = 1;
    #pragma omp parallel
    {
        #pragma omp single
        t = omp_get_num_threads();
    }
    return t;
}

static inline double deterministic_x(int seed, int j) {
    return sin(0.001 * (seed + 1) * (j + 1)) * 1000.0;
}

// ---------------- MatrixMarket reader (rank 0) ----------------

struct MmHeader {
    int rows = 0, cols = 0;
    int64_t nnz = 0;
    bool symmetric = false;
};

static inline bool starts_with(const string& s, const string& p) {
    return s.size() >= p.size() && s.compare(0, p.size(), p) == 0;
}

static MmHeader read_matrix_market_rank0(const string& path, vector<Triplet>& out) {
    ifstream file(path);
    if (!file.is_open()) throw runtime_error("Cannot open MatrixMarket file: " + path);

    MmHeader h;
    string line;
    bool symmetric = false;

    while (getline(file, line)) {
        if (line.empty()) continue;
        if (starts_with(line, "%%MatrixMarket")) {
            istringstream hs(line);
            string a,b,c,d,e;
            hs >> a >> b >> c >> d >> e;
            if (e == "symmetric") symmetric = true;
            continue;
        }
        if (line[0] == '%') continue;
        break;
    }

    require(!line.empty(), "Invalid MatrixMarket (missing size line): " + path);

    {
        istringstream ss(line);
        ss >> h.rows >> h.cols >> h.nnz;
        require((bool)ss, "Invalid size line in MatrixMarket: " + path);
    }
    h.symmetric = symmetric;

    out.clear();
    out.reserve(static_cast<size_t>(h.nnz) * (symmetric ? 2 : 1));

    int r, c;
    double v;
    for (int64_t k = 0; k < h.nnz; ++k) {
        if (!(file >> r >> c >> v)) throw runtime_error("Error reading entries: " + path);
        Triplet t{r - 1, c - 1, v};
        out.push_back(t);
        if (symmetric && t.row != t.col) out.push_back(Triplet{t.col, t.row, t.val});
    }

    if (symmetric) h.nnz = static_cast<int64_t>(out.size());
    return h;
}

// ---------------- Weak scaling synthetic generation (rank 0) ----------------

static vector<Triplet> generate_synthetic_rank0(int n, int nnz_per_row, int seed) {
    vector<Triplet> t;
    t.reserve(static_cast<size_t>(n) * static_cast<size_t>(nnz_per_row));

    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < nnz_per_row; ++k) {
            uint32_t x = static_cast<uint32_t>(i) * 2654435761u
                       ^ static_cast<uint32_t>(k) * 2246822519u
                       ^ static_cast<uint32_t>(seed) * 3266489917u;
            int j = static_cast<int>(x % static_cast<uint32_t>(n));
            double v = sin(0.01 * (seed + 1) * (i + 1) * (k + 1));
            t.push_back(Triplet{i, j, v});
        }
    }
    return t;
}

// ---------------- Distribution by owner(row)=row%P ----------------

static vector<Triplet> distribute_coo_row_cyclic(
    int rank, int P,
    const vector<Triplet>& triplets_rank0
) {
    vector<int> send_counts(P, 0), displs;
    vector<Triplet> sendbuf;

    if (rank == 0) {
        for (const auto& t : triplets_rank0) send_counts[owner_of(t.row, P)]++;

        displs.assign(P, 0);
        for (int p = 1; p < P; ++p) displs[p] = displs[p - 1] + send_counts[p - 1];

        sendbuf.resize(triplets_rank0.size());
        vector<int> cursor = displs;
        for (const auto& t : triplets_rank0) {
            int o = owner_of(t.row, P);
            sendbuf[static_cast<size_t>(cursor[o]++)] = t;
        }
    }

    int recv_count = 0;
    MPI_Scatter(send_counts.data(), 1, MPI_INT, &recv_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<Triplet> local(static_cast<size_t>(recv_count));

    MPI_Datatype MPI_TRIPLET;
    {
        int bl[3] = {1,1,1};
        MPI_Aint disp[3];
        Triplet tmp{};
        MPI_Aint base;
        MPI_Get_address(&tmp, &base);
        MPI_Get_address(&tmp.row, &disp[0]);
        MPI_Get_address(&tmp.col, &disp[1]);
        MPI_Get_address(&tmp.val, &disp[2]);
        disp[0] -= base; disp[1] -= base; disp[2] -= base;
        MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_DOUBLE};
        MPI_Type_create_struct(3, bl, disp, types, &MPI_TRIPLET);
        MPI_Type_commit(&MPI_TRIPLET);
    }

    MPI_Scatterv(
        rank == 0 ? sendbuf.data() : nullptr,
        rank == 0 ? send_counts.data() : nullptr,
        rank == 0 ? displs.data() : nullptr,
        MPI_TRIPLET,
        local.data(),
        recv_count,
        MPI_TRIPLET,
        0,
        MPI_COMM_WORLD
    );

    MPI_Type_free(&MPI_TRIPLET);

    for (auto& t : local) {
        int gr = t.row;
        require(owner_of(gr, P) == rank, "Distribution mismatch (wrong owner)");
        t.row = local_index_mod(gr, rank, P);
    }

    return local;
}

// ---------------- COO -> CSR ----------------

static bool cmp_triplet(const Triplet& a, const Triplet& b) {
    if (a.row != b.row) return a.row < b.row;
    return a.col < b.col;
}

static CsrMatrix coo_to_csr(int local_rows, int global_cols, vector<Triplet> coo_local) {
    sort(coo_local.begin(), coo_local.end(), cmp_triplet);

    CsrMatrix A;
    A.rows = local_rows;
    A.cols = global_cols;
    A.nnz  = static_cast<int64_t>(coo_local.size());

    A.row_ptr.assign(static_cast<size_t>(local_rows + 1), 0);
    A.col_ind.resize(static_cast<size_t>(A.nnz));
    A.values.resize(static_cast<size_t>(A.nnz));

    for (const auto& e : coo_local) {
        require(e.row >= 0 && e.row < local_rows, "COO invalid local row");
        A.row_ptr[static_cast<size_t>(e.row) + 1]++;
    }
    for (int i = 0; i < local_rows; ++i)
        A.row_ptr[static_cast<size_t>(i) + 1] += A.row_ptr[static_cast<size_t>(i)];

    vector<int> next = A.row_ptr;
    for (const auto& e : coo_local) {
        int pos = next[static_cast<size_t>(e.row)]++;
        A.col_ind[static_cast<size_t>(pos)] = e.col;
        A.values[static_cast<size_t>(pos)]  = e.val;
    }
    return A;
}

// ---------------- Ghost plan + remap ----------------

struct GhostPlan {
    vector<int> req_cols_packed;
    vector<int> req_pos_packed;
    vector<int> send_counts, recv_counts;
    vector<int> sdispls, rdispls;
    vector<int> recv_req_cols;
    int total_send = 0;
    int total_recv = 0;
    int n_ghosts = 0;
    int64_t bytes_indices_once = 0;
};

static GhostPlan build_ghost_plan_and_remap(
    int rank, int P,
    CsrMatrix& A,
    vector<double>& ghost_vals
) {
    GhostPlan gp;

    vector<int> ghost_cols;
    ghost_cols.reserve(static_cast<size_t>(A.nnz));
    for (int64_t k = 0; k < A.nnz; ++k) {
        int c = A.col_ind[static_cast<size_t>(k)];
        if (owner_of(c, P) != rank) ghost_cols.push_back(c);
    }
    sort(ghost_cols.begin(), ghost_cols.end());
    ghost_cols.erase(unique(ghost_cols.begin(), ghost_cols.end()), ghost_cols.end());

    gp.n_ghosts = static_cast<int>(ghost_cols.size());
    ghost_vals.assign(static_cast<size_t>(gp.n_ghosts), 0.0);

    gp.send_counts.assign(P, 0);
    for (int c : ghost_cols) gp.send_counts[owner_of(c, P)]++;
    gp.send_counts[rank] = 0;

    gp.recv_counts.assign(P, 0);
    MPI_Alltoall(gp.send_counts.data(), 1, MPI_INT, gp.recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    gp.sdispls.assign(P, 0);
    gp.rdispls.assign(P, 0);
    for (int p = 1; p < P; ++p) {
        gp.sdispls[p] = gp.sdispls[p - 1] + gp.send_counts[p - 1];
        gp.rdispls[p] = gp.rdispls[p - 1] + gp.recv_counts[p - 1];
    }

    gp.total_send = gp.sdispls[P - 1] + gp.send_counts[P - 1];
    gp.total_recv = gp.rdispls[P - 1] + gp.recv_counts[P - 1];

    gp.req_cols_packed.resize(static_cast<size_t>(gp.total_send));
    gp.req_pos_packed.resize(static_cast<size_t>(gp.total_send));

    vector<int> cursor = gp.sdispls;
    for (int pos = 0; pos < gp.n_ghosts; ++pos) {
        int c = ghost_cols[static_cast<size_t>(pos)];
        int o = owner_of(c, P);
        if (o == rank) continue;
        int idx = cursor[o]++;
        gp.req_cols_packed[static_cast<size_t>(idx)] = c;
        gp.req_pos_packed[static_cast<size_t>(idx)]  = pos;
    }

    gp.recv_req_cols.resize(static_cast<size_t>(gp.total_recv));
    MPI_Alltoallv(
        gp.req_cols_packed.data(), gp.send_counts.data(), gp.sdispls.data(), MPI_INT,
        gp.recv_req_cols.data(),   gp.recv_counts.data(), gp.rdispls.data(), MPI_INT,
        MPI_COMM_WORLD
    );

    gp.bytes_indices_once =
        static_cast<int64_t>(gp.req_cols_packed.size() + gp.recv_req_cols.size()) * (int64_t)sizeof(int);

    g_col_map.resize(static_cast<size_t>(A.nnz));
    for (int64_t k = 0; k < A.nnz; ++k) {
        int c = A.col_ind[static_cast<size_t>(k)];
        int o = owner_of(c, P);
        if (o == rank) {
            int loc = local_index_mod(c, rank, P);
            require(loc >= 0, "Owned local index negative");
            g_col_map[static_cast<size_t>(k)] = loc;
        } else {
            auto it = lower_bound(ghost_cols.begin(), ghost_cols.end(), c);
            require(it != ghost_cols.end() && *it == c, "Ghost column not found in ghost_cols");
            int pos = static_cast<int>(it - ghost_cols.begin());
            g_col_map[static_cast<size_t>(k)] = -pos - 1;
        }
    }
    return gp;
}

static inline void exchange_ghost_values(
    int rank, int P,
    const GhostPlan& gp,
    const vector<double>& x_local,
    vector<double>& ghost_vals
) {
    vector<double> send_vals(static_cast<size_t>(gp.total_recv));

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < gp.total_recv; ++i) {
        int col = gp.recv_req_cols[static_cast<size_t>(i)];
        require(owner_of(col, P) == rank, "Received request for non-owned col");
        int loc = local_index_mod(col, rank, P);
        send_vals[static_cast<size_t>(i)] = x_local[static_cast<size_t>(loc)];
    }

    vector<double> recv_vals(static_cast<size_t>(gp.total_send));

    MPI_Alltoallv(
        send_vals.data(), gp.recv_counts.data(), gp.rdispls.data(), MPI_DOUBLE,
        recv_vals.data(), gp.send_counts.data(), gp.sdispls.data(), MPI_DOUBLE,
        MPI_COMM_WORLD
    );

    for (int i = 0; i < gp.total_send; ++i) {
        int pos = gp.req_pos_packed[static_cast<size_t>(i)];
        ghost_vals[static_cast<size_t>(pos)] = recv_vals[static_cast<size_t>(i)];
    }
}

// ---------------- SpMV local ----------------

static inline void spmv_local_seq(
    const CsrMatrix& A,
    const vector<double>& x_local,
    const vector<double>& ghost_vals,
    vector<double>& y
) {
    y.assign(static_cast<size_t>(A.rows), 0.0);

    const int* rp = A.row_ptr.data();
    const double* av = A.values.data();

    for (int i = 0; i < A.rows; ++i) {
        double sum = 0.0;
        int a = rp[i], b = rp[i + 1];
        for (int k = a; k < b; ++k) {
            int m = g_col_map[static_cast<size_t>(k)];
            double xv = (m >= 0) ? x_local[static_cast<size_t>(m)]
                                 : ghost_vals[static_cast<size_t>(-m - 1)];
            sum += av[k] * xv;
        }
        y[static_cast<size_t>(i)] = sum;
    }
}

static inline void spmv_local_omp(
    const CsrMatrix& A,
    const vector<double>& x_local,
    const vector<double>& ghost_vals,
    vector<double>& y
) {
    y.assign(static_cast<size_t>(A.rows), 0.0);

    const int* rp = A.row_ptr.data();
    const double* av = A.values.data();

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < A.rows; ++i) {
        double sum = 0.0;
        int a = rp[i], b = rp[i + 1];
        for (int k = a; k < b; ++k) {
            int m = g_col_map[static_cast<size_t>(k)];
            double xv = (m >= 0) ? x_local[static_cast<size_t>(m)]
                                 : ghost_vals[static_cast<size_t>(-m - 1)];
            sum += av[k] * xv;
        }
        y[static_cast<size_t>(i)] = sum;
    }
}

// ---------------- CLI ----------------

struct Args {
    string case_name = "strong";
    string matrix_path;
    string out_csv = "results/out.csv";
    int iters = 30;
    int warmup = 3;
    int seed = 1;
    int gen_n = 0;
    int nnz_per_row = 0;
    bool use_openmp = true;
    int runs = 1;
};

static Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        string s = argv[i];
        auto need = [&](const char* key) {
            require(i + 1 < argc, string("Missing value for ") + key);
            return string(argv[++i]);
        };

        if (s == "--case") a.case_name = need("--case");
        else if (s == "--matrix") a.matrix_path = need("--matrix");
        else if (s == "--out") a.out_csv = need("--out");
        else if (s == "--iters") a.iters = stoi(need("--iters"));
        else if (s == "--warmup") a.warmup = stoi(need("--warmup"));
        else if (s == "--seed") a.seed = stoi(need("--seed"));
        else if (s == "--gen-n") a.gen_n = stoi(need("--gen-n"));
        else if (s == "--nnz-per-row") a.nnz_per_row = stoi(need("--nnz-per-row"));
        else if (s == "--no-omp") a.use_openmp = false;
        else if (s == "--runs") a.runs = stoi(need("--runs"));
        else throw runtime_error("Unknown argument: " + s);
    }

    require(a.case_name == "strong" || a.case_name == "weak", "Invalid --case (strong|weak)");
    if (a.case_name == "strong") require(!a.matrix_path.empty(), "Strong requires --matrix <file.mtx>");
    else {
        require(a.gen_n > 0, "Weak requires --gen-n <n>");
        require(a.nnz_per_row > 0, "Weak requires --nnz-per-row <k>");
    }
    require(a.iters > 0, "iters must be > 0");
    require(a.warmup >= 0, "warmup must be >= 0");
    require(a.runs > 0, "runs must be > 0");
    return a;
}

// ---------------- CSV append ----------------

static void append_csv(
    const string& out_path,
    const string& case_name,
    const string& matrix_name,
    int P, int threads,
    int run_id,
    int M, int N, int64_t nnz,
    int iters, int warmup,
    double t_total_ms, double t_comm_ms, double t_comp_ms,
    double gflops,
    int nnz_min, int nnz_max, double nnz_avg,
    double comm_bytes_per_iter
) {
    bool write_header = false;
    {
        ifstream chk(out_path);
        write_header = (!chk.good() || chk.peek() == ifstream::traits_type::eof());
    }

    ofstream out(out_path, ios::app);
    require(out.is_open(), "Cannot open CSV: " + out_path);

    if (write_header) {
        out << "case,matrix,P,threads,run_id,M,N,nnz,iters,warmup,"
               "t_total_ms,t_comm_ms,t_comp_ms,gflops,"
               "nnz_min,nnz_max,nnz_avg,comm_bytes_per_iter\n";
    }

    out << case_name << ","
        << matrix_name << ","
        << P << ","
        << threads << ","
        << run_id << ","
        << M << ","
        << N << ","
        << nnz << ","
        << iters << ","
        << warmup << ","
        << fixed << setprecision(6)
        << t_total_ms << ","
        << t_comm_ms << ","
        << t_comp_ms << ","
        << gflops << ","
        << nnz_min << ","
        << nnz_max << ","
        << nnz_avg << ","
        << comm_bytes_per_iter
        << "\n";
}

// ---------------- main ----------------

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0, P = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    try {
        Args args = parse_args(argc, argv);

        int M = 0, N = 0;
        int64_t nnz = 0;

        vector<Triplet> triplets_rank0;
        string matrix_name;

        if (rank == 0) {
            if (args.case_name == "strong") {
                MmHeader h = read_matrix_market_rank0(args.matrix_path, triplets_rank0);
                M = h.rows; N = h.cols; nnz = h.nnz;
                matrix_name = basename_noext(args.matrix_path);
            } else {
                M = args.gen_n;
                N = args.gen_n;
                triplets_rank0 = generate_synthetic_rank0(args.gen_n, args.nnz_per_row, args.seed);
                nnz = static_cast<int64_t>(triplets_rank0.size());
                matrix_name = "random";
            }
        }

        MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&nnz, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

        int local_rows = local_count_mod(M, rank, P);
        int local_cols = local_count_mod(N, rank, P);

        vector<Triplet> coo_local = distribute_coo_row_cyclic(rank, P, triplets_rank0);
        if (rank == 0) vector<Triplet>().swap(triplets_rank0);

        CsrMatrix A = coo_to_csr(local_rows, N, std::move(coo_local));

        int nnz_loc = static_cast<int>(A.nnz);
        int nnz_min = 0, nnz_max = 0;
        long long nnz_sum = 0;
        MPI_Reduce(&nnz_loc, &nnz_min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&nnz_loc, &nnz_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        long long nnz_loc_ll = nnz_loc;
        MPI_Reduce(&nnz_loc_ll, &nnz_sum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        vector<double> x_local(static_cast<size_t>(local_cols), 0.0);
        for (int jloc = 0; jloc < local_cols; ++jloc) {
            int jg = rank + jloc * P;
            x_local[static_cast<size_t>(jloc)] = deterministic_x(args.seed, jg);
        }

        vector<double> ghost_vals;
        GhostPlan gp = build_ghost_plan_and_remap(rank, P, A, ghost_vals);

        vector<double> y_local;

        auto do_spmv = [&]() {
            if (args.use_openmp) spmv_local_omp(A, x_local, ghost_vals, y_local);
            else spmv_local_seq(A, x_local, ghost_vals, y_local);
        };

        double comm_bytes_values_per_iter =
            static_cast<double>(gp.total_send + gp.total_recv) * sizeof(double);
        double comm_bytes_indices_amort =
            (args.iters > 0) ? (static_cast<double>(gp.bytes_indices_once) / args.iters) : 0.0;
        double comm_bytes_per_iter = comm_bytes_values_per_iter + comm_bytes_indices_amort;

        for (int run = 0; run < args.runs; ++run) {
            for (int w = 0; w < args.warmup; ++w) {
                MPI_Barrier(MPI_COMM_WORLD);
                exchange_ghost_values(rank, P, gp, x_local, ghost_vals);
                do_spmv();
            }

            double t_comm_sum = 0.0, t_comp_sum = 0.0, t_total_sum = 0.0;

            for (int it = 0; it < args.iters; ++it) {
                MPI_Barrier(MPI_COMM_WORLD);
                double t0 = MPI_Wtime();

                double c0 = MPI_Wtime();
                exchange_ghost_values(rank, P, gp, x_local, ghost_vals);
                double c1 = MPI_Wtime();

                double p0 = MPI_Wtime();
                do_spmv();
                double p1 = MPI_Wtime();

                double t1 = MPI_Wtime();

                t_comm_sum  += (c1 - c0);
                t_comp_sum  += (p1 - p0);
                t_total_sum += (t1 - t0);
            }

            double t_comm_max = 0.0, t_comp_max = 0.0, t_total_max = 0.0;
            MPI_Reduce(&t_comm_sum,  &t_comm_max,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&t_comp_sum,  &t_comp_max,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&t_total_sum, &t_total_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

            if (rank == 0) {
                double iters_d = static_cast<double>(args.iters);
                double t_comm_ms  = (t_comm_max  / iters_d) * 1000.0;
                double t_comp_ms  = (t_comp_max  / iters_d) * 1000.0;
                double t_total_ms = (t_total_max / iters_d) * 1000.0;

                double flops = 2.0 * static_cast<double>(nnz);
                double sec = t_total_ms / 1000.0;
                double gflops = (sec > 0.0) ? (flops / sec / 1e9) : 0.0;

                int threads = omp_threads_active();
                double nnz_avg = static_cast<double>(nnz_sum) / static_cast<double>(P);

                append_csv(args.out_csv,
                           args.case_name, matrix_name,
                           P, threads,
                           run + 1,
                           M, N, nnz,
                           args.iters, args.warmup,
                           t_total_ms, t_comm_ms, t_comp_ms,
                           gflops,
                           nnz_min, nnz_max, nnz_avg,
                           comm_bytes_per_iter);

                cout << "DONE "
                     << "run=" << (run + 1) << "/" << args.runs
                     << " case=" << args.case_name
                     << " matrix=" << matrix_name
                     << " P=" << P
                     << " T=" << threads
                     << " total_ms=" << fixed << setprecision(3) << t_total_ms
                     << " comm_ms=" << t_comm_ms
                     << " comp_ms=" << t_comp_ms
                     << " gflops=" << gflops
                     << "\n";
            }
        }

        MPI_Finalize();
        return 0;

    } catch (const exception& e) {
        if (rank == 0) cerr << "Error: " << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}

