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
// Created by Adrián on 7/4/22.
//
#include "triple_bwt.hpp"


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

bool sortby_pso(const spo_triple &a, const spo_triple &b)
{
    if(get<1>(a) == get<1>(b)){
        if(get<0>(a) == get<0>(b)){
            return get<2>(a) < get<2>(b);
        }
        return get<0>(a) < get<0>(b);
    }
    return get<1>(a) < get<1>(b);
}


bool sortby_interval(const pair<uint64_t, uint64_t> &a,
                     const pair<uint64_t, uint64_t> &b)
{
    return (a.second - a.first < b.second - b.first);
}

int main(int argc, char **argv)
{

    std::string data_file = argv[1];
    std::string so_file = argv[2];

    std::string out_data = data_file + ".sort-so";
    std::string out_so   = so_file + ".sort-so";

    vector<spo_triple> D;
    vector<std::pair<uint64_t, uint64_t>> intervals;
    unordered_map<uint64_t, uint64_t> so_hashtable;
    std::ifstream data_ifs(data_file);
    uint64_t s, p , o;
    do {
        data_ifs >> s >> p >> o;
        D.emplace_back(spo_triple(s, p, o));
    } while (!data_ifs.eof());
    data_ifs.close();
    std::sort(D.begin(), D.end(), sortby_pso);

    uint64_t p1, s1;
    uint64_t lb, line = 0;
    p1 = s1 = lb = 0;
    while(line < D.size()){
        s = get<0>(D[line]);
        p = get<1>(D[line]);
        o = get<2>(D[line]);
        if(line > 0 && (p != p1 || s != s1)){
            intervals.emplace_back(lb, line-1);
            lb = line;
        }
        p1 = p;
        s1 = s;
        ++line;
    }
    if (line> 0) intervals.emplace_back(lb, line-1);
    /*do {
        data_ifs >> s >> p >> o;
        if(data_ifs.eof()){
            if (line> 0) intervals.emplace_back(lb, line-1);
            break;
        }
        if((line > 0 && (p != p1 || o != o1))){
            intervals.emplace_back(lb, line-1);
            lb = line;
        }
        D.emplace_back(spo_triple(s, p, o));
        p1 = p;
        o1 = o;
        ++line;
    } while (true);*/


    std::cout << "D:         " << D.size() << std::endl;
    std::cout << "Intervals: " << intervals.size() << std::endl;

    std::sort(intervals.begin(), intervals.end(), sortby_interval);

    std::cout << "Intervals are sorted." << std::endl;
    uint64_t so = 1;
    for(const auto& interval : intervals){
        for(auto i = interval.first; i <= interval.second; ++i){
            o = get<2>(D[i]);
            if(so_hashtable.find(o) == so_hashtable.end()){
                so_hashtable.insert({o, so});
                ++so;
            }
        }
    }

    std::cout << "Check subjects." << std::endl;
    for(const auto& interval : intervals){
        s = get<0>(D[interval.first]);
        if(so_hashtable.find(s) == so_hashtable.end()){
            so_hashtable.insert({s, so});
            ++so;
        }
    }
    intervals.clear();
    intervals.shrink_to_fit();

    for(auto& d : D){
        s = so_hashtable[get<0>(d)];
        o = so_hashtable[get<2>(d)];
        get<0>(d) = s;
        get<2>(d) = o;
    }
    std::cout << "Triples are updated." << std::endl;

    std::sort(D.begin(), D.end(), sortby_triple);

    std::cout << "Triples are sorted." << std::endl;

    std::ofstream data_ofs(out_data);
    for(auto& d : D){
        data_ofs << get<0>(d) << " " << get<1>(d) << " " << get<2>(d) << std::endl;
    }
    data_ofs.close();

    std::cout << "Triples are printed." << std::endl;
    D.clear();
    D.shrink_to_fit();

    std::ifstream so_ifs(so_file);
    std::ofstream so_ofs(out_so);
    uint64_t id;
    std::string url, data;
    while (std::getline(so_ifs, data)) {
        auto space = data.find(' ');
        id = std::stoull(data.substr(0, space));
        url = data.substr(space+1);
        so_ofs << so_hashtable[id] << " " << url  << std::endl;
    }
    so_ifs.close();
    so_ofs.close();

    std::cout << "Dictionary is printed." << std::endl;




}