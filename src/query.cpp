//
// Created by Adrián on 2/5/23.
//

extern "C" {
    #include <matrix.h>
}
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

#define N 958844164
#define S 5420 // 1 to 5419
#define V 296008192 // 1 to...


std::string file_name(const uint i, const std::string &index){
    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << i;
    return index + "/" + ss.str() + ".mat";
}

int main(int argc, char **argv) {

    if (argc < 3) {
        std::cerr << "  Usage: " << argv[0] << " <dataset> <queries-file>" << std::endl;
        exit(1);
    }

    std::string dataset = argv[1];
    std::string queries = argv[2];
    std::string index   = dataset + ".matrices";

    //Load matrices
    matrix M[S];
    uint64_t space = 0;
    uint i;
    FILE *f;
    printf ("Reading 5420 matrices..."); fflush(stdout);
    for (i=1;i<S;i++) {
        std::string file = file_name(i, index);
        f = fopen(file.c_str(), "r");
        M[i] = matLoad(f);
        fclose(f);
        space += matSpace(M[i]);
    }
    printf (" done... %li total words (%0.2f bpt)\n",space,space*(w/8)/(float)N);

    /*std::ifstream ifs_SO(dataset + ".SO", std::ifstream::in);
    std::ifstream ifs_P(dataset + ".P", std::ifstream::in);
    std::ifstream ifs_q(argv[2], std::ifstream::in);
    //std::ofstream ofs_SO("SO", std::ofstream::out);

    std::unordered_map<string, uint64_t> map_SO;
    std::unordered_map<string, uint64_t> map_P;

    uint64_t id;
    string s_aux, data;

    while (std::getline(ifs_SO, data)) {
        auto space = data.find(' ');
        id = std::stoull(data.substr(0, space));
        s_aux = data.substr(space + 1);
        map_SO[s_aux] = id;
    }

    while (std::getline(ifs_P, data)) {
        auto space = data.find(' ');
        id = std::stoull(data.substr(0, space));
        s_aux = data.substr(space + 1);
        map_P[s_aux] = id;
    }


    typedef struct {
        uint64_t id;
        uint64_t beg;
        uint64_t end;
    } rpq_predicate;*/

}