//
// Created by Adrián on 2/5/23.
//


#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <iomanip>
#include <chrono>
#include <array>
#include "rpq_solver.hpp"
extern "C"{
#include "utilstime.h"
}


#define N 958844164
#define S 5420 // 1 to 5419
#define V 296008192 // 1 to...




bool parse_query(std::string &line,
                 std::unordered_map<std::string, uint64_t> &map_SO,
                 std::unordered_map<std::string, uint64_t> map_P,
                 std::string &query,
                 bool &flag_s,
                 int &s_id,
                 bool &flag_o,
                 int &o_id){

    std::string s_aux;
    std::stringstream X(line);
    if (line.at(0) == '?') {
        flag_s = false;
        X >> s_aux;
        if (line.at(line.size() - 3) == '?') {
            flag_o = false;
            line.pop_back();
            line.pop_back();
            line.pop_back();
        } else {
            flag_o = true;
            std::string s_aux_2;
            uint64_t i = line.size() - 1;

            while (line.at(i) != ' ') i--;
            i++;
            while (i < line.size() - 1/*line.at(i) != '>'*/)
                s_aux_2 += line.at(i++);

            if (map_SO.find(s_aux_2) != map_SO.end()) {
                o_id = map_SO[s_aux_2];
                i = 0;
                //ofs_SO << o_id << " " << s_aux_2 << endl;
                while (i < s_aux_2.size() + 1) {
                    line.pop_back();
                    i++;
                }
            } else {
                return false;
            }
        }
    } else {
        flag_s = true;
        X >> s_aux;
        if (map_SO.find(s_aux) != map_SO.end()) {
            s_id = map_SO[s_aux];
            //ofs_SO << s_id << " " << s_aux << endl;
        } else {
            return false;
        }
        if (line.at(line.size() - 3) == '?') {
            flag_o = false;
            line.pop_back();
            line.pop_back();
            line.pop_back();
        } else {
            flag_o = true;
            std::string s_aux_2;
            uint64_t i = line.size() - 2;
            while (line.at(i) != ' ') i--;
            i++;
            while (i < line.size() - 1 /*line.at(i) != '>'*/)
                s_aux_2 += line.at(i++);

            if (map_SO.find(s_aux_2) != map_SO.end()) {
                o_id = map_SO[s_aux_2];
                //ofs_SO << o_id << " " << s_aux_2 << endl;
                i = 0;
                while (i < s_aux_2.size() + 1) {
                    line.pop_back();
                    i++;
                }
            } else {
                return false;
            }
        }
    }

    std::stringstream X_a(line);
    X_a >> s_aux;
    X_a >> s_aux;
    do {
        for (uint64_t i = 0; i < s_aux.size(); i++) {
            if (s_aux.at(i) == '<') {
                std::string s_aux_2, s_aux_3;
                s_aux_2.clear();
                while (s_aux.at(i) != '>') {
                    s_aux_2 += s_aux.at(i++);
                }
                s_aux_2 += s_aux.at(i);
                if (s_aux_2[1] == '%') {
                    //Not supported
                    //return false;
                    s_aux_3 = "<" + s_aux_2.substr(2, s_aux_2.size() - 1);
                } else {
                    s_aux_3 = s_aux_2;
                }

                if (map_P.find(s_aux_3) != map_P.end()) {
                    //std::cout << s_aux_3 << " id: " << map_P[s_aux_3] << std::endl;
                    query += s_aux_2;
                } else {
                    //cout << q << ";0;0" << endl;
                    return false;
                }
            } else {
                if (s_aux.at(i) != '/' and s_aux.at(i) != ' ') {
                    if (s_aux.at(i) == '^')
                        query += '%';
                    else {
                        query += s_aux.at(i);
                    }
                }
            }
        }
    } while (X_a >> s_aux);
    return true;
}

std::string remove_unnecessary_parentheses(const std::string &query){
    std::array<uint64_t, 128> open;
    std::array<bool, 128> remove;
    uint64_t pos = 0;
    std::vector<uint64_t> pos_to_remove;
    uint64_t i = 0;
    while(i < query.size()){
        if(query[i] == '('){
            open[pos] = i;
            remove[pos] = true;
            ++pos;
        }else if (query[i] == ')'){
            --pos;
            if(i+1 < query.size() && (query[i+1] == '+' || query[i+1] == '?' || query[i+1] == '*')){
                remove[pos] = false;
                ++i;
            }
            if(remove[pos]){
                pos_to_remove.push_back(open[pos]);
                pos_to_remove.push_back(i);
            }
        }else if (query[i] == '|'){
            if(pos > 0) remove[pos-1] = false;
        }
        ++i;
    }
    if(pos_to_remove.empty()) return query;
    std::sort(pos_to_remove.begin(), pos_to_remove.end());
    pos = 0;
    std::string res;
    for(i = 0; i < query.size(); ++i){
        if(pos < pos_to_remove.size() && pos_to_remove[pos] == i){
            ++pos;
        }else{
            res += query[i];
        }
    }
    return res;

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
    rpq::solver solver(dataset, index);

    //std::string q = "?x <http://www.wikidata.org/prop/direct/P937>/<http://www.wikidata.org/prop/direct/P625> ?y ";
    //std::string q = "?x <http://www.wikidata.org/prop/direct/P680>|<http://www.wikidata.org/prop/direct/P681>|<http://www.wikidata.org/prop/direct/P682> ?y ";
    //std::string q = "?x <http://www.wikidata.org/prop/direct/P937>/<%http://www.wikidata.org/prop/direct/P625> ?y ";
    //std::string q = "?x <http://www.wikidata.org/prop/direct/P31>/(<http://www.wikidata.org/prop/direct/P279>)* ?y ";
    //std::string q = "?x (<http://www.wikidata.org/prop/direct/P194>/<http://www.wikidata.org/prop/direct/P31>)/(<http://www.wikidata.org/prop/direct/P279>)* ?y ";
    //std::string q = "?x ((<http://www.wikidata.org/prop/direct/P31>|<http://www.wikidata.org/prop/direct/P279>)/(<http://www.wikidata.org/prop/direct/P31>|<http://www.wikidata.org/prop/direct/P279>)) ?y ";
    //std::string q = "?x (<http://www.wikidata.org/prop/direct/P31>|<http://www.wikidata.org/prop/direct/P279>) ?y ";
    //std::string q = "?x (<http://www.wikidata.org/prop/direct/P31>|<http://www.wikidata.org/prop/direct/P279>)/<%http://www.wikidata.org/prop/direct/P31> ?y ";
    //std::string q = "?x <http://www.wikidata.org/prop/direct/P1448>|<http://www.wikidata.org/prop/direct/P1581>|<http://www.wikidata.org/prop/direct/P1621>|<http://www.wikidata.org/prop/direct/P242>|<http://www.wikidata.org/prop/direct/P15>|<http://www.wikidata.org/prop/direct/P18> ?y ";
    //std::string q = "?x (((<http://www.wikidata.org/prop/direct/P131>/<http://www.wikidata.org/prop/direct/P131>)/<http://www.wikidata.org/prop/direct/P131>)/<http://www.wikidata.org/prop/direct/P131>)/<http://www.wikidata.org/prop/direct/P131> ?y ";
    std::string query, line, l2;
    uint64_t i = 1;
    uint64_t elems;
    std::ifstream ifs_q(queries);
    do{
        getline(ifs_q, line);
        l2 = line;
        query.clear();

        bool flag_s, flag_o;
        int s_id, o_id;
        bool ok = parse_query(line, solver.map_SO, solver.map_P, query, flag_s, s_id, flag_o, o_id);
        if(!ok){
            std::cout << i << ";0;0" << std::endl;
        }else{
            query = remove_unnecessary_parentheses(query);
            std::cerr << query << std::endl;

            user_beg();
            // auto t1 = std::chrono::high_resolution_clock::now();
            typename rpq::solver::data_type res;
            bool rem = false;
            if (!flag_o && !flag_s) {
                res = solver.solve_var_to_var(query, rem);
            } else if (flag_o && !flag_s) {
                res = solver.solve_var_to_con(query, o_id, rem);
            } else if (!flag_o && flag_s) {
                res = solver.solve_con_to_var(query, s_id, rem);
            } else{
                //m = solver.solve_con_to_con(query, s_id, o_id, rem);
                std::cout << i << ";" << 0 << ";" << 0 << std::endl;
                ++i;
                continue;
            }
            if(res.is_transposed){
                s_matrix m = matTranspose(res.m);
                user_end();
                std::cout << i << ";" << m.elems << ";" << user_diff() << std::endl;
            }else{
                user_end();
                std::cout << i << ";" << res.m->elems << ";" << user_diff() << std::endl;
                if(res.is_tmp) matDestroy(res.m);
            }
        }
        ++i;
    }while(!ifs_q.eof());
    ifs_q.close();
}