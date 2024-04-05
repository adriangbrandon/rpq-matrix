//
// Created by Adrián on 3/4/24.
//

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include <cmath>

extern "C" {
#include <matrix.h>
}

#define SIZE 5420

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


    std::vector<matrix> m_matrices(SIZE);
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

    std::unordered_map<int, uint> table;
    uint64_t sum = 0;
    for(i = 1; i <SIZE; ++i){
        double dens = m_matrices[i]->elems/ (double) (m_matrices[i]->width * m_matrices[i]->height);
        int exp = std::ceil(log10(dens));
        auto it = table.find(exp);
        if(it == table.end()){
            table.insert({exp, 1});
        }else{
            ++it->second;
        }
        sum += m_matrices[i]->elems;
    }

    std::cout << "Avg: " << sum/(double) SIZE << std::endl;
    for(auto p : table){
        std::cout << p.first << ": " << p.second << std::endl;
    }

}