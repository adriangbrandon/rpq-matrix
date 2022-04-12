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



typedef struct {
    uint64_t p;
    uint64_t index;
    uint64_t count;

} p_data;

bool sortby_pdata(const p_data &a, const p_data &b)
{
    return a.count < b.count;
}


void BFS_node(const uint64_t p, const uint64_t s,
              std::unordered_map<uint64_t, uint64_t> &hash_table,
              std::vector<std::vector<uint64_t>> &subjects){

    auto it = hash_table.find(p);
    if(it == hash_table.end()){
        subjects.emplace_back();
        subjects.back().push_back(s);
        hash_table.insert({p, subjects.size()-1});
    }else{
        subjects[it->second].push_back(s);
    }
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

bool sortby_pos(const spo_triple &a, const spo_triple &b)
{
    if(get<1>(a) == get<1>(b)){
        if(get<2>(a) == get<2>(b)){
            return get<0>(a) < get<0>(b);
        }
        return get<2>(a) < get<2>(b);
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

    std::string out_data = data_file + ".sort-pos-bfs";
    std::string out_so   = so_file + ".sort-pos-bfs";


    std::cout << "Reading file... " << std::flush;
    vector<spo_triple> D;
    unordered_map<uint64_t, uint64_t> so_hashtable;
    std::ifstream data_ifs(data_file);
    uint64_t s, p , o;
    do {
        data_ifs >> s >> p >> o;
        D.emplace_back(spo_triple(s, p, o));
    } while (!data_ifs.eof());
    data_ifs.close();
    //std::sort(D.begin(), D.end(), sortby_pos);

    std::cout << "[done]" << std::endl;
    std::cout << "Building adjacency lists... " << std::flush;
    typedef vector<std::pair<uint64_t, uint64_t>> adj_list_type;
    vector<adj_list_type> adj_lists;
    std::unordered_map<uint64_t, uint64_t> hash_table_adj_lists;
    std::vector<uint64_t> vec_s;

    uint64_t s1, lb, line, max_so;
    line = s1 = lb = max_so = 0;
    while(line < D.size()){
        s = get<0>(D[line]);
        p = get<1>(D[line]);
        o = get<2>(D[line]);
        if(s > max_so) max_so = s;
        if(o > max_so) max_so = o;
        if(s != s1){
            hash_table_adj_lists.insert({s, adj_lists.size()});
            adj_lists.emplace_back();
            vec_s.push_back(s);
        }
        adj_lists.back().emplace_back(p, o);
        s1 = s;
        ++line;
    }
    std::cout << "[done]" << std::endl;
    std::cout << "BFS traversal... " << std::flush;
    std::unordered_map<uint64_t, uint64_t> hash_table;
    std::vector<std::vector<uint64_t>> subjects;

    //BFS traversal
    std::vector<bool> visited(max_so+1, false);
    std::queue<uint64_t> nodes;
    s = vec_s[0];
    p = 0;
    nodes.push(s);
    visited[s] = true;
    auto n_visited = 1;
    BFS_node(p, s, hash_table, subjects);
    while(!nodes.empty()){
        auto n = nodes.front();
        nodes.pop();
        auto it = hash_table_adj_lists.find(n);
        if(it != hash_table_adj_lists.end()){
            for(const auto& nn : adj_lists[it->second]){
                if(!visited[nn.second]){
                    nodes.push(nn.second);
                    BFS_node(nn.first, nn.second, hash_table, subjects);
                    visited[nn.second]=true;
                    ++n_visited;
                }
            }
        }

        if(nodes.empty()){
            std::cout << "Visited: " << n_visited << " total: " << vec_s.size() << std::endl;
            auto i = 0;
            while(i < vec_s.size() && visited[vec_s[i]]){
                ++i;
            }
            if(i >= vec_s.size()) break;
            s = vec_s[i];
            p = 0;
            nodes.push(s);
            BFS_node(p, s, hash_table, subjects);
            visited[s] = true;
            ++n_visited;
        }
    }
    vec_s.clear();
    std::cout << "[done]" << std::endl;
    std::cout << "Sorting by predicate... " << std::flush;
    std::vector<p_data> data_vec;
    data_vec.reserve(hash_table.size());
    for(const auto &it : hash_table){
        data_vec.push_back({it.first, it.second, subjects[it.second].size()});
    }
    hash_table.clear();

    std::sort(data_vec.begin(), data_vec.end(), sortby_pdata);
    std::cout << "[done]" << std::endl;
    std::cout << "Creating BFS ids... " << std::flush;
    uint64_t bfs_id = 1;
    for(const auto &d : data_vec){
        for(auto i = 0; i < d.count; ++i){
            hash_table.insert({subjects[d.index][i], bfs_id});
            ++bfs_id;
        }
    }
    data_vec.clear();
    subjects.clear();
    std::cout << "[done]" << std::endl;

    std::cout << "Updating triples... " << std::flush;
    for(auto& d : D){
        s = hash_table[get<0>(d)];
        o = hash_table[get<2>(d)];
        get<0>(d) = s;
        get<2>(d) = o;
    }
    std::cout << "[done]" << std::endl;
    std::cout << "Sorting triples... " << std::flush;

    std::sort(D.begin(), D.end(), sortby_triple);

    std::cout << "[done]" << std::endl;
    std::cout << "Writing triples... " << std::flush;
    std::ofstream data_ofs(out_data);
    for(auto& d : D){
        data_ofs << get<0>(d) << " " << get<1>(d) << " " << get<2>(d) << std::endl;
    }
    data_ofs.close();

    D.clear();
    D.shrink_to_fit();
    std::cout << "[done]" << std::endl;

    std::cout << "Writing dictionary... " << std::flush;
    std::ifstream so_ifs(so_file);
    std::ofstream so_ofs(out_so);
    uint64_t id;
    std::string url;
    do{
        so_ifs >> id >> url;
        so_ofs << hash_table[id] << " " << url  << std::endl;
    }while(!so_ifs.eof());
    so_ifs.close();
    so_ofs.close();

    std::cout << "[done]" << std::endl;




}