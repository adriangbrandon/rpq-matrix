//
// Created by Adrián on 29/3/24.
//

#include <iostream>
#include <chrono>
#include <vector>
extern "C"{
#include "matrix.h"
#include "utilstime.h"
}

int main(int argc, char **argv) {

    std::string base_name = argv[1];

    std::vector<matrix> m_matrices(20);
    FILE *f;
    uint64_t space = 0;
    printf ("Reading matrices..."); fflush(stdout);
    for (uint i=1;i<= 20 ;++i) {
        std::string file_name = base_name + "." + std::to_string(i) + ".mat";
        std::cout << file_name << std::endl;
        f = fopen(file_name.c_str(), "r");
        m_matrices[i-1] = matLoad(f);
        fclose(f);
        space += matSpace(m_matrices[i-1]);
    }
    printf (" done\n");

    matrix M;

    uint64_t sum = 0;

    std::cout << "Sum..." << std::flush;
    auto t1 = std::chrono::high_resolution_clock::now();
    for(uint i = 1; i < 20; ++i){
        time_beg();
        M = matSum(m_matrices[i-1], m_matrices[i]);
        sum += M->elems;
        matDestroy(M);
    }
    auto t2 =  std::chrono::high_resolution_clock::now();
    auto sum_ns = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
    double avg_sum = sum_ns / 19.0;
    std::cout << " done. [" << sum << "]" << std::endl;

    std::cout << "Sum RC..." << std::flush;
    sum = 0;
    t1 = std::chrono::high_resolution_clock::now();
    for(uint i = 1; i < 20; ++i){
        time_beg();
        M = matSum1(fullSide, m_matrices[i-1], m_matrices[i], fullSide);
        sum += M->elems;
        matDestroy(M);
    }
    t2 =  std::chrono::high_resolution_clock::now();
    auto sum_rec_ns = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
    double avg_sum_rec = sum_rec_ns / 19.0;
    std::cout << " done. [" << sum << "]" << std::endl;

    std::cout << "And..." << std::flush;
    sum = 0;
    t1 = std::chrono::high_resolution_clock::now();
    for(uint i = 1; i < 20; ++i){
        time_beg();
        M = matAnd(m_matrices[i-1], m_matrices[i]);
        sum += M->elems;
        matDestroy(M);
    }
    t2 =  std::chrono::high_resolution_clock::now();
    auto intersect_ns = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
    double avg_intersect = intersect_ns/19.0;
    std::cout << " done. [" << sum << "]" << std::endl;

    std::cout << "Mult..." << std::flush;
    sum = 0;
    t1 = std::chrono::high_resolution_clock::now();
    for(uint i = 1; i < 20; ++i){
        time_beg();
        M = matMult(m_matrices[i-1], m_matrices[i]);
        sum += M->elems;
        matDestroy(M);
    }
    t2 =  std::chrono::high_resolution_clock::now();
    auto product_ns = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
    double avg_product = product_ns/19.0;
    std::cout << " done. [" << sum << "]" << std::endl;

    std::cout << "Closure..." << std::flush;
    sum = 0;
    t1 = std::chrono::high_resolution_clock::now();
    for(uint i = 0; i < 20; ++i){
        time_beg();
        M = matClos(m_matrices[i], 1);
        sum += M->elems;
        matDestroy(M);
    }
    t2 =  std::chrono::high_resolution_clock::now();
    auto closure_ns = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
    double avg_closure = closure_ns / 20.0;
    std::cout << " done. [" << sum << "]" << std::endl;

    std::cerr << avg_sum << ";" << avg_sum_rec << ";" << avg_intersect << ";"<< avg_product << ";" << avg_closure << std::endl;

}