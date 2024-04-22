//
// Created by Adrián on 2/5/23.
//

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <iomanip>
#include <chrono>
#include <array>
#include <utils.hpp>

#include <bm_k2_tree_p.hpp>

#include <rpq_solver.hpp>


#define N 958844164
#define S 5420 // 1 to 5419
#define V 296008192 // 1 to...

int main(int argc, char **argv) {

    if (argc < 5) {
        std::cerr << "\tUsage: " << argv[0] << " <dataset> <queries> <n_preds> <n_triples>" << std::endl;
        exit(1);
    }

    std::string dataset = argv[1];
    std::string queries = argv[2];
    std::string index   = dataset + ".matrices";
    uint n_preds = atoi(argv[3]);
    uint n_triples = atoi(argv[4]);
    rpq::run_query<bm_k2_tree_p::wrapper>(dataset, index, queries, n_preds, n_triples);

}