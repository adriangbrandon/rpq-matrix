//
// Created by Adrián on 29/3/24.
//
#include <iomanip>

#include <iostream>
#include <chrono>
#include <vector>
extern "C"{
#include "matrix.h"
#include "utilstime.h"
}

#define NANO_TO_MILLI 1000000.0

uint64_t inline pos_matrix(uint64_t i, uint64_t d){
    return (d-1)*20+(i-1);
}

int main(int argc, char **argv) {

    std::string folder = argv[1];
    std::string type_density = argv[2];
    std::cout << folder << std::endl;

    std::vector<matrix> m_matrices(80);
    FILE *f;
    uint64_t space = 0;
    printf ("Reading matrices..."); fflush(stdout);

    uint side = 1000;
    uint size = side * side;
    uint div = 10;

    for(uint d = 1; d <= 4; ++d){
        uint elems = size / div;
        for (uint i=1;i<= 20 ;++i) {
            std::string file_name = folder + "uniform." + std::to_string(elems) + "." + std::to_string(i) + ".mat";
            std::cout << file_name << std::endl;
            f = fopen(file_name.c_str(), "r");
            uint64_t pos = pos_matrix(i, d);
            m_matrices[pos] = matLoad(f);
            fclose(f);
            space += matSpace(m_matrices[pos]);
        }
        div = div * 10;
    }

    printf (" done\n");

    matrix M;
    uint64_t sum = 0;

    //Matrix 10^-1 OP Matrix 10^-2, 10^-3, 10^-4
    uint64_t left_matrix;
    std::vector<uint> densities;

    if(type_density =="dense"){
        left_matrix = pos_matrix(1, 1);
        densities = {2, 3, 4};
    }else{
        left_matrix = pos_matrix(1, 4);
        densities = {1, 2, 3};
    }

    std::cerr << std::setprecision(8);

    for(auto d : densities){
        std::cout << "Sum..." << std::flush;
        auto t1 = std::chrono::high_resolution_clock::now();
        for(uint i = 1; i <= 20; ++i){
            M = matSum(m_matrices[left_matrix], m_matrices[pos_matrix(i, d)]);
            sum += M->elems;
            matDestroy(M);
        }
        auto t2 =  std::chrono::high_resolution_clock::now();
        auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
        std::cerr << ns / (double) (20*NANO_TO_MILLI) << std::endl;
        std::cout << " done. [" << sum << "]" << std::endl;

        std::cout << "Sum RC..." << std::flush;
        sum = 0;
        t1 = std::chrono::high_resolution_clock::now();
        for(uint i = 1; i <= 20; ++i){
            M = matOr(m_matrices[left_matrix], m_matrices[pos_matrix(i, d)]);
            sum += M->elems;
            matDestroy(M);
        }
        t2 =  std::chrono::high_resolution_clock::now();
        ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
        std::cerr << ns / (double) (20*NANO_TO_MILLI) << std::endl;
        std::cout << " done. [" << sum << "]" << std::endl;

        std::cout << "And..." << std::flush;
        sum = 0;
        t1 = std::chrono::high_resolution_clock::now();
        for(uint i = 1; i <= 20; ++i){
            M = matAnd(m_matrices[left_matrix], m_matrices[pos_matrix(i, d)]);
            sum += M->elems;
            matDestroy(M);
        }
        t2 =  std::chrono::high_resolution_clock::now();
        ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
        std::cerr << ns / (double) (20*NANO_TO_MILLI) << std::endl;
        std::cout << " done. [" << sum << "]" << std::endl;

        std::cout << "Mult..." << std::flush;
        sum = 0;
        t1 = std::chrono::high_resolution_clock::now();
        for(uint i = 1; i < 20; ++i){
            M = matMult(m_matrices[left_matrix], m_matrices[pos_matrix(i, d)]);
            sum += M->elems;
            matDestroy(M);
        }
        t2 =  std::chrono::high_resolution_clock::now();
        ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
        std::cerr << ns / (double) (20*NANO_TO_MILLI) << std::endl;
        std::cout << " done. [" << sum << "]" << std::endl;

    }



}