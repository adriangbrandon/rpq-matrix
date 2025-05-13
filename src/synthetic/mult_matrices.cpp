//
// Created by adrian on 8/5/25.
//
#include <iostream>
#include <chrono>
#include <vector>
#include <bm_k2_tree.hpp>

#define NANO_TO_MILLI 1000000.0

int main(int argc, char **argv) {

    std::string base_name = argv[1];

    typedef bm_k2_tree::wrapper wrapper;
    typedef typename wrapper::matrix_type matrix;

    std::vector<matrix> m_matrices(20);
    FILE *f;
    uint64_t space = 0;
    printf ("Reading matrices..."); fflush(stdout);

    for (uint i=0; i<m_matrices.size(); i++) {
        std::string file = base_name + "." + std::to_string(i+1) + ".mat";
        f = fopen(file.c_str(), "r");
        m_matrices[i] = wrapper::load(f);
        fclose(f);
        space += wrapper::space(m_matrices[i]);
    }


    matrix tmp;
    uint64_t sum = 0;

    std::cout << "Mult..." << std::flush;
    sum = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    for(uint i = 1; i < 20; ++i){
        tmp = wrapper::mult(m_matrices[i-1], m_matrices[i]);
        sum += tmp->elems;
    }
    auto t2 =  std::chrono::high_resolution_clock::now();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
    std::cerr << ns / (double) (19*NANO_TO_MILLI) << std::endl;
    std::cout << " done. [" << sum << "]" << std::endl;


}