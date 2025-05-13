//
// Created by adrian on 8/5/25.
//

//
// Created by adrian on 8/5/25.
//
#include <iostream>
#include <chrono>
#include <vector>
#include <bm_k2_tree.hpp>
#include <file.hpp>
#include <algorithm>

#define NANO_TO_MILLI 1000000.0

typedef std::vector<std::vector<bool>> matrix_type;

void print_points(std::vector<std::pair<uint32_t, uint32_t>> &vp) {
    std::sort(vp.begin(), vp.end());
    for (auto & p : vp) {
        std::cout << "Point: (" << p.first << ", " << p.second << ")" << std::endl;
    }
}

void summ (matrix_type &r, matrix_type &a, matrix_type &b){
    auto n = a.size();
    auto m = a[0].size(); //    auto p = b[0].size();
    r.resize(n);
    for(auto & row : r) {
        row.resize(m);
        //row.resize(p);
    }
    for(auto i = 0; i < n; ++i) {
        for(auto j = 0; j < m; ++j) {
            r[i][j] = a[i][j] || b[i][j];
        }
    }
}

void points2matrix(std::vector<std::pair<uint32_t, uint32_t>> &points, matrix_type &m) {
    for(auto &p : points) {
        m[p.first][p.second] = true;
    }
}



int main(int argc, char **argv) {

    std::string matrices_dir = argv[1];

    auto files = util::file::read_directory(matrices_dir);
    std::sort(files.begin(), files.end());

    typedef bm_k2_tree::wrapper wrapper;
    typedef typename wrapper::matrix_type matrix;

    std::vector<matrix> m_matrices(files.size());
    uint64_t space = 0;
    std::cout << "Reading matrices..." << std::flush;

    FILE *f;
    uint nmatrices = 0;
    for(auto & file : files) {
        std::string file_path = matrices_dir + "/" + file;
        std::cout << "[" << file_path << "]" << std::endl;
        f = fopen(file_path.c_str(), "r");
        m_matrices[nmatrices] = wrapper::load(f);
        fclose(f);
        space += wrapper::space(m_matrices[nmatrices]);
        ++nmatrices;
    }
    std::cout << " done. [" << space << " B]" << std::endl;

    matrix tmp;
    uint64_t sum = 0;

    
    std::cout << "Sum..." << std::flush;
    sum = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    for(uint i = 1; i < nmatrices; ++i){
        tmp = wrapper::sum(m_matrices[i-1], m_matrices[i]);
        sum += tmp->elems;
    }
    auto t2 =  std::chrono::high_resolution_clock::now();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
    std::cerr << space << ";" << ns / (double) ((nmatrices-1)*NANO_TO_MILLI) << std::endl;
    std::cout << " done. [" << sum << "]" << std::endl;


}
