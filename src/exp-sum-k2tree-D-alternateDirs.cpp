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
#include <cmath>
#include <iomanip>

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

    if (argc <4) {
        std::cerr << "Usage: " << argv[0] << " <matrices_dir1> <matrices_dir1> <max-matrices-to-read>" << std::endl;
        return 0;
    }
    std::string matrices_dir1 = argv[1];
    std::string matrices_dir2 = argv[2];
    uint32_t max_matrices_to_read = atoi(argv[3]);

    //  Reads up to max_matrices_to_read matrices from 1st directory
    auto files1 = util::file::read_directory(matrices_dir1);
    std::sort(files1.begin(), files1.end());

    //  Reads up to max_matrices_to_read matrices from 1st directory
    auto files2 = util::file::read_directory(matrices_dir2);
    std::sort(files2.begin(), files2.end());

    uint32_t nmatrices = std::min(files1.size(), files2.size());
    nmatrices = std::min(max_matrices_to_read, nmatrices);
    std::cout << nmatrices << " operations :: M1_i · M2_i will be performed" << std::endl;

    uint32_t i; ////////////
    typedef bm_k2_tree::wrapper wrapper;
    typedef typename wrapper::matrix_type matrix;

    std::vector<matrix> m_matrices1(files1.size());
    uint64_t space = 0,aux;
    std::cout << "Reading matrices M1_i from 1st dir ..." << std::endl<<std::flush;

    FILE *f1;
    i= 0;
    for(auto & file : files1) {
        std::string file_path = matrices_dir1 + "/" + file;
        //std::cout << "[" << file_path << "]" << std::endl;
        f1 = fopen(file_path.c_str(), "r");
        m_matrices1[i] = wrapper::load(f1);
        fclose(f1);
        aux = wrapper::space(m_matrices1[i]) *8;  // 64bit words to bytes
        std::cout << "\t[M1_" << std::left << std::setw((int)std::log10(nmatrices-1) + 1 )  << i << std::setfill(' ')
                  << ": " << file_path << "] bytes: " << aux << std::endl;
        space += aux;
        ++i;
        if (i>=nmatrices) break;
    }

    std::vector<matrix> m_matrices2(files2.size());
    std::cout << "Reading matrices M2_i from 2nd dir ..." << std::endl<<std::flush;

    FILE *f2;
    i= 0;
    bool special_case=false;
    for(auto & file : files2) {
        if ((nmatrices == 1) && (files2.size()>1) && !special_case) {
            //SPECIAL CASE: only 1 operation will be RUN, and we are operating matrices within the same dir
            //Then, if possible -ie. there are more than 1 matrix in dir2- we simply skip the first matrix from dir2
            //so that the second matrix is set at m_matrices2[0].
            // This permits to avoid M·M (operating the first matrix in dir1=dir2 by itself).
            // and operates the first matrix (in dir1) by the second matrix (in dir2=dir2) instead. ;)
            special_case=true;
            continue;
        }

        std::string file_path = matrices_dir2 + "/" + file;
        f2 = fopen(file_path.c_str(), "r");
        //NOTE: If dir1 and dir2 are the same directory
        //      and there are at least 2 matrices (nmatrices >1)
        //      we make M1[i] · M2[(i+1)%nmatrices] to avoid M·M (operating the same matrix by itself).
        //      Otherwise (only 1 matrix in directory) we will operating M·M :(
        m_matrices2[(i+1)%nmatrices] = wrapper::load(f2);
          
        fclose(f2);
        aux= wrapper::space(m_matrices2[(i+1)%nmatrices]) *8; // 64bit words to bytes
        std::cout << "\t[M2_" << std::left << std::setw((int)std::log10(nmatrices-1) + 1 )  << (i+1)%nmatrices << std::setfill(' ')
                  << ": " << file_path << "] bytes: " << aux << std::endl;
        space += aux;
        ++i;
        if (i>=nmatrices) break;
    }
    //space *=8;   //from words to bytes

    std::cout << " Loads completed: [" << space << " Bytes (sum of all "<< nmatrices*2<< " matrices)]" << std::endl;

    matrix tmp;
    uint64_t sum = 0;

    std::cout << "Sum..." << std::flush;
    sum = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    for(uint i = 0; i <nmatrices; ++i){
        tmp = wrapper::sum(m_matrices1[i], m_matrices2[i]);
        sum += tmp->elems;
    }
    auto t2 =  std::chrono::high_resolution_clock::now();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
    std::cerr << space/(nmatrices*2) << ";" << ns / (double) ((nmatrices)*NANO_TO_MILLI) << std::endl;
    std::cout << " done. [ results: " << sum << "] " << "; ["<<nmatrices<< " calls to sum]; (avg-space per matrix used= " << space/(nmatrices*2) <<")" << std::endl;

}
