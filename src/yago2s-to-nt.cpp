/***
BSD 2-Clause License

Copyright (c) 2018, Adrián
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/


//
// Created by Adrián on 13/6/22.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include "triple_bwt.hpp"


int main(int argc, char **argv){

    if(argc != 3){
        std::cout << argv[0] << " <input> <output>";
        return 0;
    }

    std::string input  = argv[1];
    std::string output = argv[2];

    std::ifstream ifs_SO(output + ".SO", std::ifstream::in);
    std::ifstream ifs_P(output + ".P", std::ifstream::in);
    //std::ofstream ofs_SO("SO", std::ofstream::out);

    std::unordered_map<uint64_t, string> map_SO;
    std::unordered_map<uint64_t, string> map_P;

    uint64_t id;
    string s_aux, data;

    while (std::getline(ifs_SO, data)) {
        auto space = data.find(' ');
        id = std::stoull(data.substr(0, space));
        s_aux = data.substr(space+1);
        map_SO[id] = s_aux;
    }

    while (std::getline(ifs_P, data)) {
        auto space = data.find(' ');
        id = std::stoull(data.substr(0, space));
        s_aux = data.substr(space+1);
        map_P[id] = s_aux;
    }

    uint64_t s, p, o;
    std::ifstream in(input);
    std::ofstream out(output);
    while(true){
        in >> s >> p >> o;
        if(in.eof()) break;
        out << map_SO[s] << " " << map_P[p] << " " << map_SO[o] << " ." << std::endl;
    }
    in.close();
    out.close();

    return 1;

}