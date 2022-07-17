#include <iostream>
#include "triple_bwt_v3.hpp"
#include <fstream>
#include <sdsl/construct.hpp>

using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;


int main(int argc, char **argv) {

    if (argc < 3) {
        cerr << "  Usage: " << argv[0] << " <ring-index-file> <queries-file>" << endl;
        exit(1);
    }

    ring_rpq_bfs graph;

    graph.load(string(argv[1]));

    cout << "  Index loaded " << (float) graph.size() / (graph.n_triples() / 2) << " bytes per triple" << endl;

    std::ifstream ifs_SO(string(argv[1]) + ".SO", std::ifstream::in);
    std::ifstream ifs_P(string(argv[1]) + ".P", std::ifstream::in);
    std::ifstream ifs_q(argv[2], std::ifstream::in);
    //std::ofstream ofs_SO("SO", std::ofstream::out);

    std::unordered_map<string, uint64_t> map_SO;
    std::unordered_map<string, uint64_t> map_P;

    uint64_t id;
    string s_aux, data;

    while (std::getline(ifs_SO, data)) {
        auto space = data.find(' ');
        id = std::stoull(data.substr(0, space));
        s_aux = data.substr(space+1);
        map_SO[s_aux] = id;
    }

    while (std::getline(ifs_P, data)) {
        auto space = data.find(' ');
        id = std::stoull(data.substr(0, space));
        s_aux = data.substr(space+1);
        map_P[s_aux] = id;
    }


    typedef struct {
        uint64_t id;
        uint64_t beg;
        uint64_t end;
    } rpq_predicate;

    std::unordered_map<std::string, uint64_t> pred_map;
    std::string query, line;
    uint64_t s_id, o_id, n_line = 0, q = 0;
    bool flag_s, flag_o, skip_flag;
    std::vector<std::pair<uint64_t, uint64_t>> query_output;
    std::vector<word_t> B_array(4 * graph.n_labels(), 0);
    std::vector<word_t> B_array2(4 * graph.n_labels(), 0);

    high_resolution_clock::time_point start, stop;
    double total_time = 0.0;
    duration<double> time_span;

    uint64_t n_predicates, n_operators;
    bool is_negated_pred, is_a_path, is_or;
    std::string query_type;
    uint64_t first_pred_id, last_pred_id;

    uint64_t bound = 1000000; // bound for the number of results

    do {
        getline(ifs_q, line);
        pred_map.clear();
        query.clear();
        query_type.clear();
        query_output.clear();
        flag_s = flag_o = false;
        skip_flag = false;
        n_line++;

        stringstream X(line);
        q++;

        n_predicates = 0;
        is_negated_pred = false;
        n_operators = 0;
        is_a_path = true;
        is_or = true;

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
                string s_aux_2;
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
                } else
                    skip_flag = true;

            }
        } else {
            flag_s = true;
            X >> s_aux;
            if (map_SO.find(s_aux) != map_SO.end()) {
                s_id = map_SO[s_aux];
                //ofs_SO << s_id << " " << s_aux << endl;
            } else
                skip_flag = true;

            if (line.at(line.size() - 3) == '?') {
                flag_o = false;
                line.pop_back();
                line.pop_back();
                line.pop_back();
            } else {
                flag_o = true;
                string s_aux_2;
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
                } else
                    skip_flag = true;
            }
        }

        stringstream X_a(line);
        X_a >> s_aux;

        if (!skip_flag) {
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
                            s_aux_3 = "<" + s_aux_2.substr(2, s_aux_2.size() - 1);
                            is_negated_pred = true;
                        } else {
                            s_aux_3 = s_aux_2;
                        }

                        if (map_P.find(s_aux_3) != map_P.end()) {
                            query += s_aux_2;
                            last_pred_id = pred_map[s_aux_3] = map_P[s_aux_3];
                            n_predicates++;
                            if (n_predicates == 1)
                                first_pred_id = pred_map[s_aux_3];
                        } else {
                            cout << q << ";0;0" << endl;
                            skip_flag = true;
                            break;
                        }
                    } else {
                        if (s_aux.at(i) != '/' and s_aux.at(i) != ' ') {
                            n_operators++;
                            if (s_aux.at(i) == '^')
                                query += '%';
                            else {
                                query += s_aux.at(i);
                                if (s_aux.at(i) != '(' and s_aux.at(i) != ')') query_type += s_aux.at(i);
                                is_a_path = false;
                                if (s_aux.at(i) != '|') is_or = false;
                            }
                        } else if (s_aux.at(i) == '/') {
                            n_operators++;
                            query_type += s_aux.at(i);
                        }
                    }
                }
            } while (!skip_flag and X_a >> s_aux);

            if (!skip_flag) {
                start = high_resolution_clock::now();
                if (!flag_s and !flag_o) {
                    if (query_type.size() == 0)
                        graph.rpq_var_to_var_so(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                                n_operators, false);
                    else if (query_type == "/*")
                        /*graph.rpq_var_to_var_os(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                                n_operators, false);*/
                        graph.rpq_var_to_var_so(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                                n_operators, false);
                    else if (query_type == "+" or query_type == "*")
                        graph.rpq_var_to_var_so(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                                n_operators, false);
                    else if (query_type == "|")
                        graph.or_query_var_to_var(query, 2, bound, pred_map, query_output);
                    else if (query_type == "||")
                        graph.or_query_var_to_var(query, 3, bound, pred_map, query_output);
                    else if (is_or and n_predicates > 3)
                        graph.or_query_var_to_var(query, n_predicates, bound, pred_map, query_output);
                    else if (query_type[query_type.size() - 1] == '*' and query_type[0] != '*')
                        /*graph.rpq_var_to_var_os(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                                n_operators, false);*/
                        graph.rpq_var_to_var_so(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                                n_operators, false);
                    else if (query_type[0] == '*' or query_type[0] == '+')
                        /*graph.rpq_var_to_var_so(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                                n_operators, false);*/
                        graph.rpq_var_to_var_os(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                                n_operators, false);
                    else if (graph.pred_selectivity(first_pred_id) <= graph.pred_selectivity(last_pred_id))
                        graph.rpq_var_to_var_so(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                                n_operators, is_a_path);
                        /*graph.rpq_var_to_var_split(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                                   n_operators, is_a_path, bound);*/
                    else
                        graph.rpq_var_to_var_os(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                                n_operators, is_a_path);
                    /*graph.rpq_var_to_var_split(query, pred_map, B_array, query_output, n_predicates, is_negated_pred,
                                           n_operators, is_a_path, bound);*/
                } else {
                    if (flag_s) {
                        graph.rpq_const_s_to_var_o(query, pred_map, B_array, s_id, query_output, n_predicates,
                                                   is_negated_pred, n_operators, is_a_path);

                    } else {
                        graph.rpq_var_s_to_const_o(query, pred_map, B_array, o_id, query_output, n_predicates,
                                                   is_negated_pred, n_operators, is_a_path);
                    }
                }
                stop = high_resolution_clock::now();
                time_span = duration_cast<microseconds>(stop - start);
                total_time = time_span.count();

                cout << q << ";" << query_output.size() << ";" << (uint64_t) (total_time * 1000000000ULL) << endl;
            } else skip_flag = false;
        } else {
            cout << q << ";0;0" << endl;
            skip_flag = false;
        }
    } while (!ifs_q.eof());

    return 0;
}

