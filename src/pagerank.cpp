/* PageRank for graph matrices represented as k2-trees (C++ version).
 *
 * Port of k2pagerank.c using the k2-tree operations in
 * k2_tree_operations.hpp.  The only structural operation needed is the
 * matrix-vector product z = M*y, which is provided by cds::op::mult_vec().
 *
 * As in the original program:
 *   - the matrix stored on disk is the TRANSPOSE of the adjacency matrix,
 *     with self-loops already removed, so that the right product M*y
 *     implements the PageRank left multiplication;
 *   - col_count_file holds one uint32_t per node with its out-degree
 *     (= number of nonzeros in the corresponding column of the original
 *     matrix).  Nodes with out-degree 0 are dangling nodes.
 *
 * Iteration (N nodes, damping d, current rank X):
 *   for i: if outd[i]==0: dnr += X[i]   else Y[i] = X[i]/outd[i]
 *   Z = M*Y
 *   for i: Z[i] = d*Z[i] + (d*dnr + (1-d))/N
 *   X = Z
 * We start with X = (1/N, ..., 1/N) and stop after maxiter iterations or
 * when the L1 difference between consecutive rank vectors is < eps.
 * As in the original, only Y and Z are kept; X is recovered from Y at the end.
 *
 */
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>
#include <unistd.h>            // getopt

#include "bm_k2_tree.hpp"

// ---------------------------------------------------------------------------
//  top-k via a min-heap (ported verbatim from k2pagerank.c)
// ---------------------------------------------------------------------------
static void swap_int(int *a, int *b) { int t = *a; *a = *b; *b = t; }

static void minHeapify(const double v[], int arr[], int n, int i) {
    int smallest = i, left = 2 * i + 1, right = 2 * i + 2;
    if (left  < n && v[arr[left]]  < v[arr[smallest]]) smallest = left;
    if (right < n && v[arr[right]] < v[arr[smallest]]) smallest = right;
    if (smallest != i) {
        swap_int(&arr[i], &arr[smallest]);
        minHeapify(v, arr, n, smallest);
    }
}

// fills arr[0..k-1] with the indices of the k largest values of v[0..n-1]
static void kLargest(const double v[], int arr[], int n, int k) {
    for (int i = 0; i < k; ++i) arr[i] = i;
    for (int i = k / 2 - 1; i >= 0; --i) minHeapify(v, arr, k, i);
    for (int i = k; i < n; ++i)
        if (v[i] > v[arr[0]]) { arr[0] = i; minHeapify(v, arr, k, 0); }
}

static void usage_and_exit(const char *name) {
    fprintf(stderr, "Usage:\n\t%s [options] matrix col_count_file\n", name);
    fprintf(stderr, "\t\t-v          verbose (repeat for per-iteration log)\n");
    fprintf(stderr, "\t\t-m maxiter  maximum number of iterations (default 100)\n");
    fprintf(stderr, "\t\t-e eps      stop if error<eps (default: ignore error)\n");
    fprintf(stderr, "\t\t-d df       damping factor (default 0.9)\n");
    fprintf(stderr, "\t\t-k K        show top K nodes (default 3)\n");
    exit(1);
}

int main(int argc, char **argv) {
    time_t start_wc = time(NULL);

    int maxiter = 100, topk = 3, verbose = 0;
    double dampf = 0.9, eps = -1;

    int c;
    opterr = 0;
    while ((c = getopt(argc, argv, "m:e:d:k:v")) != -1) {
        switch (c) {
            case 'v': verbose++;            break;
            case 'm': maxiter = atoi(optarg); break;
            case 'e': eps     = atof(optarg); break;
            case 'd': dampf   = atof(optarg); break;
            case 'k': topk    = atoi(optarg); break;
            default:  usage_and_exit(argv[0]);
        }
    }
    if (argc - optind != 2)             usage_and_exit(argv[0]);
    if (maxiter < 1 || topk < 1)        usage_and_exit(argv[0]);
    if (dampf < 0 || dampf > 1)         usage_and_exit(argv[0]);

    const char *matrix_file = argv[optind];
    const char *ccol_file   = argv[optind + 1];

    // ----------- read column-count file (out-degrees), get N -----------------
    std::vector<uint32_t> outd;
    {
        std::ifstream f(ccol_file, std::ios::binary | std::ios::ate);
        if (!f) { fprintf(stderr, "Cannot open col_count_file '%s'\n", ccol_file); return 1; }
        std::streamsize bytes = f.tellg();
        size_t n = (size_t)(bytes / (std::streamsize)sizeof(uint32_t));
        if (n < 1) { fprintf(stderr, "Invalid col_count_file\n"); return 1; }
        outd.resize(n);
        f.seekg(0);
        f.read(reinterpret_cast<char *>(outd.data()), (std::streamsize)(n * sizeof(uint32_t)));
        if (!f) { fprintf(stderr, "Cannot read col_count_file\n"); return 1; }
    }
    const size_t N = outd.size();

    // ----------- load the k2-tree matrix -------------------------------------

    typedef bm_k2_tree::wrapper wrapper;
    typedef typename wrapper::matrix_type matrix;

    matrix A;
    FILE* f = fopen(matrix_file, "r");
    A = wrapper::load(f);
    if (A->height < N) {
        fprintf(stderr, "Matrix size (%llu) < col_count size (%llu)\n",
                (unsigned long long) A->height, (unsigned long long)N);
        return 1;
    }

    if (verbose) {
        long dn = 0, arcs = 0;
        for (size_t i = 0; i < N; ++i)
            if (outd[i] == 0) dn++; else arcs += (long)outd[i];
        fprintf(stderr, "Number of nodes: %zu\n", N);
        fprintf(stderr, "Number of dangling nodes: %ld\n", dn);
        fprintf(stderr, "Number of arcs: %ld\n", arcs);
    }

    // ----------- init rank / aux vectors -------------------------------------
    // x_0 = (1/N, ..., 1/N); we store directly y_0 and the dangling-rank sum.
    double* y = new double[N];
    double dnr = 0;
    for (size_t i = 0; i < N; ++i) {
        if (outd[i] == 0) dnr += (y[i] = 1.0 / (double)N);
        else              y[i] = (1.0 / (double)N) / outd[i];
    }
    // ----------- main iteration ----------------------------------------------
    int iter = 0;
    double delta = 11 + eps;   // ensure we don't stop before the first iteration
    while (iter < maxiter && delta >= eps) {
        double* z;
        z = wrapper::mult_vec(A, y);                 // z = M * y  (z is (re)sized to N)

        /*for (uint64_t i = 0; i < 100; ++i) {
            std::cout << iter << " " << z[i] << std::endl;
        }
        exit(0);*/

        double teleport = (dnr * dampf + 1 - dampf) / (double)N;
        dnr = delta = 0;
        for (size_t i = 0; i < N; ++i) {
            double nextri = dampf * z[i] + teleport;
            if (outd[i] == 0) {
                delta += fabs(nextri - y[i]);          // dangling: X[i] == y[i]
                dnr   += (y[i] = nextri);
            } else {
                delta += fabs(nextri - y[i] * outd[i]); // X[i] == y[i]*outd[i]
                y[i]   = nextri / outd[i];
            }
        }
        ++iter;
        if (verbose > 1) fprintf(stderr, "Iteration %d, delta=%g\n", iter, delta);
    }

    // recover the actual rank vector x from y
    for (size_t i = 0; i < N; ++i)
        if (outd[i] != 0) y[i] = y[i] * outd[i];
    double* x = y;   // x is the rank vector now

    if (verbose) {
        if (delta > eps) fprintf(stderr, "Stopped after %d iterations, delta=%g\n", iter, delta);
        else             fprintf(stderr, "Converged after %d iterations, delta=%g\n", iter, delta);
        double sum = 0;
        for (size_t i = 0; i < N; ++i) sum += x[i];
        fprintf(stderr, "Sum of ranks: %f (should be 1)\n", sum);
    }

    // ----------- top-k nodes --------------------------------------------------
    if ((size_t)topk > N) topk = (int)N;
    std::vector<int> top(topk), aux(topk);
    kLargest(x, aux.data(), (int)N, topk);   // NOTE: assumes N < 2^31
    for (int i = topk - 1; i >= 0; --i) {           // extract in decreasing order
        top[i]  = aux[0];
        aux[0]  = aux[i];
        minHeapify(x, aux.data(), i, 0);
    }

    if (verbose) {
        fprintf(stderr, "Top %d ranks:\n", topk);
        for (int i = 0; i < topk; ++i)
            fprintf(stderr, "  %d %lf\n", top[i], x[top[i]]);
    }
    fprintf(stdout, "Top:");
    for (int i = 0; i < topk; ++i) fprintf(stdout, " %d", top[i]);
    fprintf(stdout, "\n");

    fprintf(stderr, "Elapsed time: %.0lf secs\n", (double)(time(NULL) - start_wc));
    return 0;
}
