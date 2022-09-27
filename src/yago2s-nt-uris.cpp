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

bool is_literal(const std::string &uri){
    return (uri.find("<\"") == 0);
}

std::string get_literal(const std::string &literal){
    return literal.substr(1, literal.size()-2);
}

bool is_uri(const std::string &uri){
    return (uri.find("<http://") == 0);
}

std::string get_prefix(const std::string &uri){
    auto pos = uri.find(':');
    if(pos!=std::string::npos){
        return uri.substr(1, pos);
    }else{
        return "";
    }
}

void get_id_uri(const std::string &line, std::string &id, std::string &uri){
    auto pos = line.find(" ");
    id = line.substr(0, pos);
    uri = line.substr(pos+1);
}

std::pair<bool, std::string> fix_prefix(const std::string &uri, const std::string &prefix, const std::string &p_uri){
    auto pos = uri.find(prefix);
    if(pos == std::string::npos) return {false, ""};
    return {true, std::regex_replace(uri, std::regex(prefix), p_uri)};
}

std::string base(const std::string &uri, const std::string &uri_base){
    auto str_aux = uri.substr(1, uri.size()-2);
    return "<" + uri_base + str_aux + ">";
}


std::string transform(const std::string &uri,
                      const std::vector<std::pair<std::string, std::string>>& prefixes,
                      const std::string &uri_base){

    if(is_uri(uri)){
        return uri;
    }else if (is_literal(uri)){
        return get_literal(uri);
    }else{
        for(const auto &p : prefixes){
            auto fix = fix_prefix(uri, p.first, p.second);
            if(fix.first){
                return fix.second;
            }
        }
    }
    return base(uri, uri_base);



}

int main(int argc, char **argv){

    if(argc != 3){
        std::cout << argv[0] << " <input> <output>";
        return 0;
    }

    std::string input  = argv[1];
    std::string output = argv[2];

    std::string so_f = input + ".SO";
    std::string p_f = input + ".P";

    std::string so_o_f = output + ".SO";
    std::string p_o_f = output + ".P";

    std::string uri_base = "http://yago-knowledge.org/resource/";
    std::pair<std::string, std::string> rdf = {"rdf:", "http://www.w3.org/1999/02/22-rdf-syntax-ns#"};
    std::pair<std::string, std::string> rdfs = {"rdfs:", "http://www.w3.org/2000/01/rdf-schema#"};
    std::pair<std::string, std::string> owl = {"owl:", "http://www.w3.org/2002/07/owl#"};
    std::pair<std::string, std::string> xsd = {"xsd:", "http://www.w3.org/2001/XMLSchema#"};
    std::pair<std::string, std::string> skos = {"skos:", "http://www.w3.org/2004/02/skos/core#"};

    std::vector<std::pair<std::string, std::string>> prefixes = {rdf, rdfs, owl, xsd, skos};

    std::ifstream so_in(so_f);
    std::ofstream so_out(so_o_f);
    std::string line;
    std::unordered_map<std::string, uint64_t> so_map;
    std::unordered_map<std::string, uint64_t> p_map;
    std::string id, uri;
    std::vector<spo_triple> D;
    while(std::getline(so_in, line)){
        get_id_uri(line, id, uri);
        so_out << id << " " << transform(uri, prefixes, uri_base) << std::endl;
    }
    so_in.close();
    so_out.close();

    std::ifstream p_in(p_f);
    std::ofstream p_out(p_o_f);
    while(std::getline(p_in, line)){
        get_id_uri(line, id, uri);
        p_out << id << " " << transform(uri, prefixes, uri_base) << std::endl;
    }
    p_in.close();
    p_out.close();
    return 1;

}