//
// Created by Adrián on 3/4/24.
//

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <chrono>

extern "C" {
#include <matrix.h>
}

#define SIZE 5420

std::vector<matrix> m_matrices(SIZE);

void computeDensity(uint64_t i){
    double dens = m_matrices[i]->elems/ (double) (m_matrices[i]->width * m_matrices[i]->height);
    int exp = std::ceil(log10(dens));
    std::cout << "Matriz: " << i << std::endl;
    std::cout << "Density: " << dens << std::endl;
    std::cout << "Exp: " << exp << std::endl;
}

void computeDensity(matrix i){
    double dens = i->elems/ (double) (i->width * i->height);
    int exp = std::ceil(log10(dens));
    std::cout << "Reuslt " << std::endl;
    std::cout << "Density: " << dens << std::endl;
    std::cout << "Exp: " << exp << std::endl;
}

void product(uint64_t i, uint64_t j){
    std::cout << "M[" << i << "]/M[" << j << "]" << std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    auto M = matMult(m_matrices[i], m_matrices[j]);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
    std::cout << "Time " << i << "/" << j << " :" << ns << std::endl;
    std::cout << "Elems: " << M->elems << std::endl;
    computeDensity(i);
    std::cout << std::endl;
    computeDensity(j);
    std::cout << std::endl;
    computeDensity(M);
    std::cout << std::endl;
}

std::string file_name(const uint i, const std::string &index){
    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << i;
    return index + "/" + ss.str() + ".mat";
}


int main(int argc, char **argv) {

    if (argc < 2) {
        std::cerr << "  Usage: " << argv[0] << " <dataset>" << std::endl;
        exit(1);
    }
    std::string dataset = argv[1];
    std::string index   = dataset + ".matrices";



    uint64_t space = 0;
    uint i;
    FILE *f;
    std::cout << "Reading 5420 matrices..." << std::flush;
    for (i=1;i<SIZE;i++) {
        std::string file = file_name(i, index);
        f = fopen(file.c_str(), "r");
        m_matrices[i] = matLoad(f);
        fclose(f);
        space += matSpace(m_matrices[i]);
    }
    std::cout << "[done]" << std::endl;

    product(16, 1082);
    product(221, 2083);


}