//
// Created by Adrián on 16/4/24.
//

#ifndef RPQ_MATRIX_UTILS_HPP
#define RPQ_MATRIX_UTILS_HPP

#include <vector>
#include <algorithm>
#include <rpq_solver.hpp>

namespace rpq {

    namespace utils {

        static bool parse_query(std::string &line,
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

        static std::string remove_unnecessary_parentheses(const std::string &query){
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
    }

    template<class wrapper_type>
    void run_query(const std::string &dataset, const std::string &index, const std::string &queries, uint n_preds, uint n_triples){
        typedef rpq::solver<wrapper_type> solver_type;
        typedef typename solver_type::matrix matrix;
        typedef typename solver_type::s_matrix s_matrix;

        solver_type solver(dataset, index, n_preds, n_triples);

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
            bool ok = rpq::utils::parse_query(line, solver.map_SO, solver.map_P, query, flag_s, s_id, flag_o, o_id);
            if(!ok){
                std::cout << i << ";0;0" << std::endl;
            }else{
                query = rpq::utils::remove_unnecessary_parentheses(query);
                std::cerr << query << std::endl;

                auto t1 = std::chrono::high_resolution_clock::now();
                wrapper_type::time_begin();
                // auto t1 = std::chrono::high_resolution_clock::now();
                typename solver_type::data_type res;
                bool rem = false;
                if (!flag_o && !flag_s) {
                    res = solver.solve_var_to_var(query, rem);
                } else if (flag_o && !flag_s) {
                    res = solver.solve_var_to_con(query, o_id, rem);
                } else if (!flag_o && flag_s) {
                    res = solver.solve_con_to_var(query, s_id, rem);
                } else{
                    res = solver.solve_con_to_con(query, s_id, o_id, rem);
                }
                if(res.is_transposed){
                    s_matrix m = wrapper_type::transpose(res.m);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto t = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
                    std::cout << i << ";" << m.elems << ";" <<t << std::endl;
                }else{
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto t = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
                    std::cout << i << ";" << res.m->elems << ";" << t << std::endl;
                    if(res.is_tmp) wrapper_type::destroy(res.m);
                }
            }
            ++i;
        }while(!ifs_q.eof());
        ifs_q.close();
    }
}
#endif //RPQ_MATRIX_UTILS_HPP
