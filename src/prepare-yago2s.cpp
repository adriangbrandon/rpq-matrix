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

inline bool valid_line(std::string &line){
    return !line.empty() && line[0] != '@' && line[0] != '#';
}

inline void fix_data(std::string &d){
    if(d[0] != '<'){
        d = '<' + d + '>';
    }
}

bool sortbysec(const std::pair<std::string, uint64_t> &a,
               const std::pair<std::string, uint64_t> &b)
{
    return (a.second < b.second);
}

bool sortby_triple(const spo_triple &a, const spo_triple &b)
{
    if(get<0>(a) == get<0>(b)){
        if(get<1>(a) == get<1>(b)){
            return get<2>(a) < get<2>(b);
        }
        return get<1>(a) < get<1>(b);
    }
    return get<0>(a) < get<0>(b);
}

void map_to_file(const std::string &file_name, std::unordered_map<std::string, uint64_t> &map){

    std::vector<std::pair<std::string, uint64_t>> pairs(map.begin(), map.end());
    map.clear();
    std::sort(pairs.begin(), pairs.end(), sortbysec);
    std::ofstream out(file_name);
    for(const auto &pair : pairs){
        out << pair.second << " " << pair.first << std::endl;
    }
    out.close();
}

void vector_to_file(const std::string &file_name, std::vector<spo_triple> &vec){

    std::ofstream out(file_name);
    for(const auto &spo : vec){
        out <<  get<0>(spo) << " " << get<1>(spo) << " " << get<2>(spo) << std::endl;
    }
    out.close();
}

int main(int argc, char **argv){

    if(argc != 3){
        std::cout << argv[0] << " <input.ttl> <output>";
        return 0;
    }

    std::string input  = argv[1];
    std::string output = argv[2];

    std::ifstream in(input);
    std::string line;
    std::unordered_map<std::string, uint64_t> so_map;
    std::unordered_map<std::string, uint64_t> p_map;
    uint64_t so_new_id = 1;
    uint64_t p_new_id = 1;
    std::vector<spo_triple> D;
    while (std::getline(in, line)) {
        if (valid_line(line)) {
            std::istringstream iss(line);
            std::string s, p, o, point;
            uint64_t s_id, p_id, o_id;
            std::getline( iss, s, '\t' );
            std::getline( iss, p, '\t' );
            std::getline( iss, o, '\t' );
            std::getline( iss, point, '\t' );
            //iss >> s >> p >> o >> point;
            fix_data(s);
            fix_data(p);
            fix_data(o);
            std::cout << s << " " << p << " " << o << std::endl;
            auto it_so = so_map.find(s);
            if(it_so != so_map.end()){
                s_id = it_so->second;
            }else{
                s_id = so_new_id++;
                so_map.insert({s, s_id});
            }
            auto it_p = p_map.find(p);
            if(it_p != p_map.end()){
                p_id = it_p->second;
            }else{
                p_id = p_new_id++;
                p_map.insert({p, p_id});
            }
            it_so = so_map.find(o);
            if(it_so != so_map.end()){
                o_id = it_so->second;
            }else{
                o_id = so_new_id++;
                so_map.insert({o, o_id});
            }
            D.emplace_back(spo_triple{s_id, p_id, o_id});
        }
    }
    in.close();
    std::sort(D.begin(), D.end(), sortby_triple);
    //Build file
    vector_to_file(output, D);
    //Build files .P and .SO
    map_to_file(output + ".SO", so_map);
    map_to_file(output + ".P", p_map);
    return 1;

}