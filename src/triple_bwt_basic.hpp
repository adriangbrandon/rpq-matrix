/*
 * triple_bwt_rpq.hpp
 * Copyright (C) 2021 Diego Arroyuelo
 * 
 * This is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This software is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TRIPLE_L
#define TRIPLE_L

#include <cstdint>
#include <chrono>
#include <set>
#include <unordered_set>
#include <stdio.h>
#include <stdlib.h>
#include <sdsl/init_array.hpp>

#include "bwt.hpp"
#include "bwt-C-nose.hpp"
#include "bwt_interval.hpp"
#include "utils.hpp"
#include "RpqAutomata.hpp"
#include "RpqTree.hpp"
#include "query_config.hpp"
#include "parse_query.cpp"
#include "wt_intersection.hpp"
#include "selectivity.hpp"

#define ELEMENTS 0
#define RUN_QUERY 0

#ifdef CHECK_MEM
#include <malloc_count/malloc_count.h>
#include <malloc_count/memprofile.h>
#endif

using namespace std::chrono;

typedef struct {
    bwt_interval interval;
    word_t current_D;
} interval_state_type;

typedef struct {
    uint64_t element;
    std::vector<std::pair<uint64_t, uint64_t>> solutions;
} element_solution_type;

template<class Container = std::stack<interval_state_type>>
class ring_rpq {
    bwt_nose L_S;
    bwt L_P;

    uint64_t real_max_P;

    uint64_t max_S;
    uint64_t max_P;
    uint64_t max_O;
    uint64_t nTriples;  // number of triples

public:


    ring_rpq() { ; }

    // Assumes the triples have been stored in a vector<spo_triple>
    /*ring_rpq(vector<spo_triple>& D, bool verbose = true)
    {
        uint64_t i, pos_c;
        vector<spo_triple>::iterator it, triple_begin, triple_end;
        uint64_t n;

        // for every triple, adds its reverse (using a predicate
        // shifted by max_P, so 1 becomes max_P+1 in the reverse,
        // and so on)

        uint64_t d = D.size();
        spo_triple triple_aux;

        max_S = get<0>(D[0]);
        max_P = get<1>(D[0]);
        max_O = get<2>(D[0]);

        for (i = 1; i < D.size(); i++) {
            if (max_S < get<0>(D[i])) max_S = get<0>(D[i]);
            if (max_P < get<1>(D[i])) max_P = get<1>(D[i]);
            if (max_O < get<2>(D[i])) max_O = get<2>(D[i]);
        }

        real_max_P = max_P;

        if (verbose) cout << "  > Adding inverted edges";
        for (uint64_t i = 0; i < d; i++) {
            triple_aux = D[i];
            D.push_back(spo_triple(get<2>(D[i]), get<1>(D[i]), get<0>(D[i])));
            get<0>(D[i]) = get<0>(triple_aux);
            get<1>(D[i]) = get<1>(D[i]) + real_max_P;
            get<2>(D[i]) = get<2>(triple_aux);
        }
        D.shrink_to_fit();

        if (verbose) cout << "... [done]" << endl;
        max_S = get<0>(D[0]);
        max_P = get<1>(D[0]);
        max_O = get<2>(D[0]);

        for (i = 1; i < D.size(); i++) {
            if (max_S < get<0>(D[i])) max_S = get<0>(D[i]);
            if (max_P < get<1>(D[i])) max_P = get<1>(D[i]);
            if (max_O < get<2>(D[i])) max_O = get<2>(D[i]);
        }

        {
            bit_vector bv_s(max_S + 1, 0);
            for (i = 0; i < D.size(); i++) {
                bv_s[get<0>(D[i])] = 1;
            }

            bit_vector bv_o(max_O + 1, 0);
            for (i = 0; i < D.size(); i++) {
                bv_o[get<2>(D[i])] = 1;
            }

            uint64_t _c = 0;
            for (i = 1; i < max_S + 1; i++) {
                if (!bv_s[i]) {
                    D.push_back(spo_triple(i, max_P + 1, max_O + 1));
                    _c++;
                }
            }
            //cout << _c << " nodes are no subjects" << endl;

            _c = 0;
            for (i = 1; i < max_O + 1; i++) {
                if (!bv_o[i]) {
                    D.push_back(spo_triple(max_O + 1, max_P + 1, i));
                    _c++;
                }
            }

            //cout << _c << " nodes are no objects" << endl;
        }

        max_S++;
        max_O++;
        max_P++;

        triple_begin = D.begin();
        triple_end = D.end();

        if (verbose) cout << "  > Triples set = " << D.size() * sizeof(spo_triple) << " bytes" << endl;
        fflush(stdout);


        n = nTriples = triple_end - triple_begin;

        if (verbose) cout << "  > Determining number of elements per symbol";
        fflush(stdout);
        uint64_t alphabet_SO = (max_S < max_O) ? max_O : max_S;

        std::vector<uint32_t> M_O(alphabet_SO + 1), M_P(max_P + 1);

        for (i = 0; i <= alphabet_SO; i++) M_O[i] = 0;

        for (i = 0; i <= max_P; i++) M_P[i] = 0;

        for (i = 0; i < D.size(); i++) {
            M_O[std::get<2>(D[i])]++;
            M_P[std::get<1>(D[i])]++;
        }

        M_O.shrink_to_fit();
        M_P.shrink_to_fit();

        if (verbose) cout << "... [done]\n  > Sorting out triples";
        // Sorts the triples lexycographically
        sort(triple_begin, triple_end);
        if (verbose) cout << "... [done]" << endl;

        dat = new uint32_t[3 * n + 2];
        uint32_t *t = dat;

        if (verbose) cout << "  > Generating int vector of the triples";
        fflush(stdout);
        for (i = 0, it = triple_begin; it != triple_end; it++, i++) {
            t[3 * i] = std::get<0>(*it);
            t[3 * i + 1] = std::get<1>(*it) + n;
            t[3 * i + 2] = std::get<2>(*it) + 2 * n;
        }
        t[3 * n] = 3 * n;
        t[3 * n + 1] = 0;
        D.clear();
        D.shrink_to_fit();
        if (verbose) cout << "... [done]" << endl;
        if (verbose) cout << "  > Building the suffix array";
        uint64_t *sa = new uint64_t[3 * n + 1];

        for (i = 0; i < n; i++) sa[i] = 3 * i;
        for (i = n; i < 2 * n; i++) sa[i] = 3 * (i - n) + 1;
        for (i = 2 * n; i < 3 * n; i++) sa[i] = 3 * (i - 2 * n) + 2;
        sa[3 * n] = 3 * n;
        qsort(sa + n, n, sizeof(uint64_t), compare);
        qsort(sa + 2 * n, n, sizeof(uint64_t), compare);

        if (verbose) cout << "... [done]" << endl;
        if (verbose) cout << "  > Building bwts";
        std::fstream fp;
        fp.open("bwt_s.tmp", ios::out | ios::binary);
        uint64_t tmp;

        uint64_t j;
        for (j = 1, i = n; i < 2 * n; i++) {
            if (sa[i] >= 3 * n) continue;
            if (sa[i] == 0) {
                tmp = t[3 * n - 1];
            } else
                tmp = t[sa[i] - 1];
            fp.write((char *) &tmp, sizeof(uint64_t));
            j++;
        }

        // L_p
        fp.close();
        fp.open("bwt_p.tmp", ios::out | ios::binary);
        for (j = 1, i = 2 * n; i < 3 * n; i++) {
            if (sa[i] >= 3 * n) continue;
            if (sa[i] == 0) {
                tmp = t[3 * n - 1] - n;
            } else
                tmp = t[sa[i] - 1] - n;
            fp.write((char *) &tmp, sizeof(uint64_t));
            j++;
        }
        fp.close();
        delete[] sa;
        delete[] t;
        cout << "...done" << endl;

        {
            std::unordered_map<uint64_t, uint64_t> hash_table;
            int_vector<> bwt_s(n + 1);
            bwt_s[0] = 0;
            fp.open("bwt_s.tmp", ios::in | ios::binary);
            for (j = 1; j <= n; j++) {
                fp.read((char *) &tmp, sizeof(uint64_t));
                bwt_s[j] = tmp;
                auto ht_it = hash_table.find(tmp);
                if(ht_it == hash_table.end()){
                    hash_table.insert({tmp, j});
                }else{
                    ht_it->second = j;
                }
            }
            fp.close();
            util::bit_compress(bwt_s);
            std::remove("bwt_s.tmp");

            // Then S
            uint64_t c;
            if (verbose) cout << "  > Building C_S";
            vector < uint64_t > C_S;
            uint64_t cur_pos = 1;
            C_S.push_back(0);  // Dummy value
            C_S.push_back(cur_pos);
            for (c = 2; c <= max_P; c++) {
                cur_pos += M_P[c - 1];
                C_S.push_back(cur_pos);
            }
            C_S.push_back(n + 1);
            C_S.shrink_to_fit();
            if (verbose) cout << "... [done]" << endl;
            M_P.clear();
            M_P.shrink_to_fit();
            if (verbose) {
                cout << "  > Building the wavelet tree for bwt_s";
                fflush(stdout);
            }
            L_S = bwt_nose(bwt_s, C_S);
            if (verbose) cout << "... [done]" << endl;
            C_S.clear();
            C_S.shrink_to_fit();
        }

        // Then P
        {
            int_vector<> bwt_p(n + 1);
            bwt_p[0] = 0;
            fp.open("bwt_p.tmp", ios::in | ios::binary);
            for (j = 1; j <= n; j++) {
                fp.read((char *) &tmp, sizeof(uint64_t));
                bwt_p[j] = tmp;
            }
            fp.close();
            util::bit_compress(bwt_p);
            std::remove("bwt_p.tmp");

            if (verbose) cout << "  > Building C_P";
            vector < uint64_t > C_P;
            uint64_t cur_pos = 1;
            C_P.push_back(0);  // Dummy value
            C_P.push_back(cur_pos);
            uint64_t c;
            for (c = 2; c <= alphabet_SO; c++) {
                cur_pos += M_O[c - 1];
                C_P.push_back(cur_pos);
            }
            C_P.push_back(n + 1);
            C_P.shrink_to_fit();
            if (verbose) cout << "... [done]" << endl;
            M_O.clear();
            M_O.shrink_to_fit();
            if (verbose) {
                cout << "  > Building the wavelet tree for bwt_p";
                fflush(stdout);
            }
            L_P = bwt(bwt_p, C_P);
            if (verbose) cout << "... [done]" << endl;
            C_P.clear();
            C_P.shrink_to_fit();
        }
        if (verbose) cout << "-- Index constructed successfully" << endl;
    };*/

    // Assumes the triples have been stored in a vector<spo_triple>
    ring_rpq(vector<spo_triple>& D, bool verbose = true)
    {
        uint64_t i, pos_c;
        vector<spo_triple>::iterator it, triple_begin, triple_end;
        uint64_t n;

        // for every triple, adds its reverse (using a predicate
        // shifted by max_P, so 1 becomes max_P+1 in the reverse,
        // and so on)

        uint64_t d = D.size();
        spo_triple triple_aux;

        max_S = get<0>(D[0]);
        max_P = get<1>(D[0]);
        max_O = get<2>(D[0]);

        for (i = 1; i < D.size(); i++) {
            if (max_S < get<0>(D[i])) max_S = get<0>(D[i]);
            if (max_P < get<1>(D[i])) max_P = get<1>(D[i]);
            if (max_O < get<2>(D[i])) max_O = get<2>(D[i]);
        }

        real_max_P = max_P;

        if (verbose) cout << "  > Adding inverted edges";
        for (uint64_t i = 0; i < d; i++) {
            triple_aux = D[i];
            D.push_back(spo_triple(get<2>(D[i]), get<1>(D[i]), get<0>(D[i])));
            get<0>(D[i]) = get<0>(triple_aux);
            get<1>(D[i]) = get<1>(D[i]) + real_max_P;
            get<2>(D[i]) = get<2>(triple_aux);
        }
        D.shrink_to_fit();

        if (verbose) cout << "... [done]" << endl;
        max_S = get<0>(D[0]);
        max_P = get<1>(D[0]);
        max_O = get<2>(D[0]);

        for (i = 1; i < D.size(); i++) {
            if (max_S < get<0>(D[i])) max_S = get<0>(D[i]);
            if (max_P < get<1>(D[i])) max_P = get<1>(D[i]);
            if (max_O < get<2>(D[i])) max_O = get<2>(D[i]);
        }

        {
            bit_vector bv_s(max_S + 1, 0);
            for (i = 0; i < D.size(); i++) {
                bv_s[get<0>(D[i])] = 1;
            }

            bit_vector bv_o(max_O + 1, 0);
            for (i = 0; i < D.size(); i++) {
                bv_o[get<2>(D[i])] = 1;
            }

            uint64_t _c = 0;
            for (i = 1; i < max_S + 1; i++) {
                if (!bv_s[i]) {
                    D.push_back(spo_triple(i, max_P + 1, max_O + 1));
                    _c++;
                }
            }
            //cout << _c << " nodes are no subjects" << endl;

            _c = 0;
            for (i = 1; i < max_O + 1; i++) {
                if (!bv_o[i]) {
                    D.push_back(spo_triple(max_O + 1, max_P + 1, i));
                    _c++;
                }
            }

            //cout << _c << " nodes are no objects" << endl;
        }

        max_S++;
        max_O++;
        max_P++;

        triple_begin = D.begin();
        triple_end = D.end();

        if (verbose) cout << "  > Triples set = " << D.size() * sizeof(spo_triple) << " bytes" << endl;
        fflush(stdout);


        n = nTriples = triple_end - triple_begin;

        if (verbose) cout << "  > Determining number of elements per symbol";
        fflush(stdout);
        uint64_t alphabet_SO = (max_S < max_O) ? max_O : max_S;

        std::vector<uint32_t> M_O(alphabet_SO + 1, 0), M_P(max_P + 1, 0);

        for (i = 0; i < D.size(); i++) {
            M_O[std::get<2>(D[i])]++;
            M_P[std::get<1>(D[i])]++;
        }

        M_O.shrink_to_fit();
        M_P.shrink_to_fit();

        if (verbose) cout << "... [done]\n  > Sorting out triples in OSP order";
        // Sorts the triples in OSP order
        sort(D.begin(), D.end(), [](const spo_triple& a,
                                    const spo_triple& b) { return std::tie(std::get<2>(a), std::get<0>(a), std::get<1>(a))
                                                                  < std::tie(std::get<2>(b), std::get<0>(b), std::get<1>(b));});
        if (verbose) cout << "... [done]" << endl;

        //Building BWT_P

        {
            vector<uint64_t> new_C_P;
            if (verbose) cout << "  > Building C_P";
            uint64_t cur_pos = 1;
            new_C_P.push_back(0);  // Dummy value
            new_C_P.push_back(cur_pos);
            for (uint64_t c = 2; c <= alphabet_SO; c++) {
                cur_pos += M_O[c-1];
                new_C_P.push_back(cur_pos);
            }
            new_C_P.push_back(n+1);
            new_C_P.shrink_to_fit();
            M_O.clear();
            if (verbose) cout << "... [done]" << endl;
            if (verbose) cout << "  > Building bwtp";
            int_vector<> new_P(n+1);
            new_P[0] = 0;
            for (i=1; i<=n; i++)
                new_P[i] = std::get<1>(D[i-1]);
            util::bit_compress(new_P);

            L_P = bwt(new_P, new_C_P);
            if (verbose) cout << "... [done]" << endl;
        }

        if (verbose) cout << "... [done]\n  > Stable sorting out triples in POS order";
        stable_sort(D.begin(), D.end(), [](const spo_triple& a,
                                           const spo_triple& b) {return std::get<1>(a) < std::get<1>(b);});

        //Building BWT_S
        {
            vector<uint64_t> new_C_S;
            if (verbose) cout << "  > Building C_S";
            uint64_t cur_pos = 1;
            new_C_S.push_back(0);  // Dummy value
            new_C_S.push_back(cur_pos);
            for (uint64_t c = 2; c <= max_P; c++) {
                cur_pos += M_P[c-1];
                new_C_S.push_back(cur_pos);
            }
            new_C_S.push_back(n+1);
            new_C_S.shrink_to_fit();
            M_P.clear();
            if (verbose) cout << "... [done]" << endl;
            if (verbose) cout << "  > Building bwts";
            int_vector<> new_S(n+1);
            new_S[0] = 0;
            for (i=1; i<=n; i++)
                new_S[i] = std::get<0>(D[i-1]);
            util::bit_compress(new_S);

            L_S = bwt_nose(new_S, new_C_S);
            if (verbose) cout << "... [done]" << endl;
        }


        if (verbose) cout << "-- Index constructed successfully" << endl;
    };

    // Size of the index in bytes
    uint64_t size() {
        //cout << "L_O: " << (float)L_O.size()*8/nTriples << endl;
        //cout << "L_S: " << (float)L_S.size()/nTriples << " bytes per Triple" << endl;
        //cout << "L_P: " << (float)L_P.size()/nTriples << " bytes per triple" << endl;

        //cout << "L_S.size() = " << L_S.size() << " bytes" << endl;
        //cout << "L_P.size() = " << L_P.size() << " bytes" << endl;
        //cout << "nTriples = " << nTriples << endl;
        return L_S.size() + L_P.size();
    }

    uint64_t n_triples() {
        return nTriples;
    }

    uint64_t n_labels() {
        return max_P;
    }

    void save(string filename) {
        L_S.save(filename + ".bwts");
        L_P.save(filename + ".bwtp");

        std::ofstream ofs(filename + ".nTriples");
        ofs << nTriples << endl;
        ofs << real_max_P << endl;
        ofs << max_S << endl;
        ofs << max_P << endl;
        ofs << max_O << endl;
    };

    void load(string filename) {
        //cout << "Loading L_S" << endl; fflush(stdout);
        L_S.load(filename + ".bwts");
        //cout << "Loading L_P" << endl; fflush(stdout);
        L_P.load(filename + ".bwtp");
        //cout << "Loading done" << endl; fflush(stdout);
        std::ifstream ifs(filename + ".nTriples");
        ifs >> nTriples;
        ifs >> real_max_P;
        ifs >> max_S;
        ifs >> max_P;
        ifs >> max_O;
    };

    uint64_t pred_selectivity(uint64_t pred_id) {
        return L_S.get_C(pred_id + 1) - L_S.get_C(pred_id);
    };


    inline uint64_t pred_reverse(uint64_t pred_id) const{
        return (pred_id > real_max_P) ? pred_id - real_max_P : pred_id + real_max_P;
    }



private:

    inline interval_state_type& last_element(std::stack<interval_state_type> &stack){
        return stack.top();
    }


    inline interval_state_type first_element(std::stack<interval_state_type> &stack){
        return stack.top();
    }


    inline interval_state_type& last_element(std::queue<interval_state_type> &queue){
        return queue.back();
    }

    inline interval_state_type first_element(std::queue<interval_state_type> &queue){
        return queue.front();
    }


    void next_step_objects(RpqAutomata &A, std::vector<word_t> &B_array,
                           initializable_array<word_t> &D_array,
                           word_t current_D, bwt_interval &I_p,
                           Container &ist_container,
                           std::vector<uint64_t> &objects){

        //PART1: Finding predicates from the object whose range in L_p is I_p
        std::vector<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>> pred_vec;
        L_P.all_active_p_values_in_range_test<word_t>(I_p.left(), I_p.right(), B_array, current_D, pred_vec);

        std::pair<uint64_t, uint64_t> interval;
        uint64_t c;
        //PART2: For each predicate, the values in L_S are obtained by using a backward step
        for (uint64_t i = 0; i < pred_vec.size(); i++) {
            interval = L_P.backward_step_test(I_p.left(), I_p.right(), pred_vec[i].first,
                                              pred_vec[i].second.first,
                                              pred_vec[i].second.second);
            c = L_S.get_C(pred_vec[i].first);
            std::vector<std::tuple<uint64_t, word_t, std::pair<uint64_t, uint64_t>>> subj_vec;
            L_S.all_active_s_values_in_range_test<word_t>(c + interval.first, c + interval.second,
                                                          D_array,
                                                          (word_t) A.next(current_D, pred_vec[i].first, BWD),
                                                          subj_vec);

            //PART3: Map the range of each subject to the range of objects
            for (uint64_t j = 0; j < subj_vec.size(); j++) {

                if (A.atFinal(get<1>(subj_vec[j]), BWD))
                    objects.push_back(get<0>(subj_vec[j]));

                push_merge_interval(ist_container, subj_vec[j]);
            }
        }

    }

    void next_step_objects_with_bound(RpqAutomata &A, std::vector<word_t> &B_array,
                           initializable_array<word_t> &D_array,
                           word_t current_D, bwt_interval &I_p,
                           Container &ist_container,
                           std::vector<uint64_t> &objects,
                           uint64_t bound){

        //PART1: Finding predicates from the object whose range in L_p is I_p
        std::vector<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>> pred_vec;
        L_P.all_active_p_values_in_range_test<word_t>(I_p.left(), I_p.right(), B_array, current_D, pred_vec);

        std::pair<uint64_t, uint64_t> interval;
        uint64_t c;
        //PART2: For each predicate, the values in L_S are obtained by using a backward step
        for (uint64_t i = 0; i < pred_vec.size(); i++) {
            interval = L_P.backward_step_test(I_p.left(), I_p.right(), pred_vec[i].first,
                                              pred_vec[i].second.first,
                                              pred_vec[i].second.second);
            c = L_S.get_C(pred_vec[i].first);
            std::vector<std::tuple<uint64_t, word_t, std::pair<uint64_t, uint64_t>>> subj_vec;
            L_S.all_active_s_values_in_range_test<word_t>(c + interval.first, c + interval.second,
                                                          D_array,
                                                          (word_t) A.next(current_D, pred_vec[i].first, BWD),
                                                          subj_vec);

            //PART3: Map the range of each subject to the range of objects
            for (uint64_t j = 0; j < subj_vec.size() and objects.size() < bound; j++) {

                if (A.atFinal(get<1>(subj_vec[j]), BWD))
                    objects.push_back(get<0>(subj_vec[j]));

                push_merge_interval(ist_container, subj_vec[j]);
            }
        }

    }

    void next_step_const_to_var(RpqAutomata &A, std::vector<word_t> &B_array,
                             initializable_array<word_t> &D_array,
                             word_t current_D, bwt_interval &I_p,
                             uint64_t initial_object,
                             Container &ist_container,
                             std::vector<std::pair<uint64_t, uint64_t>> &solutions,
                             bool const_to_var = true){

        //PART1: Finding predicates from the object whose range in L_p is I_p
        std::vector<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>> pred_vec;
        L_P.all_active_p_values_in_range_test<word_t>(I_p.left(), I_p.right(), B_array, current_D, pred_vec);

        std::pair<uint64_t, uint64_t> interval;
        uint64_t c;
        //PART2: For each predicate, the values in L_S are obtained by using a backward step
        for (uint64_t i = 0; i < pred_vec.size(); i++) {
            interval = L_P.backward_step_test(I_p.left(), I_p.right(), pred_vec[i].first,
                                              pred_vec[i].second.first,
                                              pred_vec[i].second.second);
            c = L_S.get_C(pred_vec[i].first);
            std::vector<std::tuple<uint64_t, word_t, std::pair<uint64_t, uint64_t>>> subj_vec;
            L_S.all_active_s_values_in_range_test<word_t>(c + interval.first, c + interval.second,
                                                          D_array, (word_t) A.next(current_D, pred_vec[i].first, BWD)
                                                          ,subj_vec);

            //PART3: Map the range of each subject to the range of objects
            for (uint64_t j = 0; j < subj_vec.size(); j++) {
                push_merge_interval(ist_container, subj_vec[j]);
                if (A.atFinal(get<1>(subj_vec[j]), BWD)) {
                    if(const_to_var){
                        solutions.emplace_back(initial_object, get<0>(subj_vec[j]));
                    }else{
                        solutions.emplace_back(get<0>(subj_vec[j]), initial_object);
                    }
                }
            }
        }

    }

    void next_step_const_to_var_verbose(RpqAutomata &A, std::vector<word_t> &B_array,
                                initializable_array<word_t> &D_array,
                                word_t current_D, bwt_interval &I_p,
                                uint64_t initial_object,
                                Container &ist_container,
                                std::vector<std::pair<uint64_t, uint64_t>> &solutions,
                                bool const_to_var = true){

        //PART1: Finding predicates from the object whose range in L_p is I_p
        std::vector<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>> pred_vec;
        //auto t0 = std::chrono::high_resolution_clock::now();
        L_P.all_active_p_values_in_range_test<word_t>(I_p.left(), I_p.right(), B_array, current_D, pred_vec);
        //auto t1 = std::chrono::high_resolution_clock::now();
        //auto part1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
        //std::cout << "1;" << part1 <<";" << pred_vec.size() << std::endl;
        std::pair<uint64_t, uint64_t> interval;
        uint64_t c;
        //PART2: For each predicate, the values in L_S are obtained by using a backward step
        for (uint64_t i = 0; i < pred_vec.size(); i++) {
            //auto t2 = std::chrono::high_resolution_clock::now();
            interval = L_P.backward_step_test(I_p.left(), I_p.right(), pred_vec[i].first,
                                              pred_vec[i].second.first,
                                              pred_vec[i].second.second);
            c = L_S.get_C(pred_vec[i].first);
            std::vector<std::tuple<uint64_t, word_t, std::pair<uint64_t, uint64_t>>> subj_vec;
            L_S.all_active_s_values_in_range_test<word_t>(c + interval.first, c + interval.second,
                                                          D_array, (word_t) A.next(current_D, pred_vec[i].first, BWD)
                    ,subj_vec);
            //auto t3 = std::chrono::high_resolution_clock::now();
            //auto part2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count();
            //std::cout << "2;" << part2 <<";" << i << ";" << interval.first << ";"
            //<< interval.second << ";" << subj_vec.size() << std::endl;
            //PART3: Map the range of each subject to the range of objects
            //t2 = std::chrono::high_resolution_clock::now();
            for (uint64_t j = 0; j < subj_vec.size(); j++) {
                push_merge_interval(ist_container, subj_vec[j]);
                if (A.atFinal(get<1>(subj_vec[j]), BWD)) {
                    if(const_to_var){
                        solutions.emplace_back(initial_object, get<0>(subj_vec[j]));
                    }else{
                        solutions.emplace_back(get<0>(subj_vec[j]), initial_object);
                    }
                }
            }
            //t3 = std::chrono::high_resolution_clock::now();
            //auto part3 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count();
            //std::cout << "3;" << part3 << std::endl;
        }

    }

    void next_step_const_to_var_with_bound(RpqAutomata &A, std::vector<word_t> &B_array,
                             initializable_array<word_t> &D_array,
                             word_t current_D, bwt_interval &I_p,
                             uint64_t initial_object,
                             Container &ist_container,
                             std::vector<std::pair<uint64_t, uint64_t>> &solutions,
                             uint64_t bound,
                             bool const_to_var = true){

        //PART1: Finding predicates from the object whose range in L_p is I_p
        std::vector<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>> pred_vec;
        L_P.all_active_p_values_in_range_test<word_t>(I_p.left(), I_p.right(), B_array, current_D, pred_vec);

        std::pair<uint64_t, uint64_t> interval;
        uint64_t c;
        //PART2: For each predicate, the values in L_S are obtained by using a backward step
        for (uint64_t i = 0; i < pred_vec.size(); i++) {
            interval = L_P.backward_step_test(I_p.left(), I_p.right(), pred_vec[i].first,
                                              pred_vec[i].second.first,
                                              pred_vec[i].second.second);
            c = L_S.get_C(pred_vec[i].first);
            std::vector<std::tuple<uint64_t, word_t, std::pair<uint64_t, uint64_t>>> subj_vec;
            L_S.all_active_s_values_in_range_test<word_t>(c + interval.first, c + interval.second,
                                                          D_array,
                                                          (word_t) A.next(current_D, pred_vec[i].first, BWD),
                                                          subj_vec);

            //PART3: Map the range of each subject to the range of objects
            for (uint64_t j = 0; j < subj_vec.size() and solutions.size() < bound; j++) {
                push_merge_interval(ist_container, subj_vec[j]);

                if (A.atFinal(get<1>(subj_vec[j]), BWD)) {
                    if (const_to_var) {
                        solutions.emplace_back(initial_object, get<0>(subj_vec[j]));
                    } else {
                        solutions.emplace_back(get<0>(subj_vec[j]), initial_object);
                    }
                }
            }
        }

    }


    bool next_step_const_to_const(RpqAutomata &A, std::vector<word_t> &B_array,
                             initializable_array<word_t> &D_array,
                             word_t current_D, bwt_interval &I_p,
                             uint64_t const_1, //start with it
                             uint64_t const_2,
                             Container &ist_container){


        //PART1: Finding predicates from the object whose range in L_p is I_p
        std::vector<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>> pred_vec;
        L_P.all_active_p_values_in_range_test<word_t>(I_p.left(), I_p.right(), B_array, current_D, pred_vec);

        std::pair<uint64_t, uint64_t> interval;
        uint64_t c;
        //PART2: For each predicate, the values in L_S are obtained by using a backward step
        for (uint64_t i = 0; i < pred_vec.size(); i++) {
            interval = L_P.backward_step_test(I_p.left(), I_p.right(), pred_vec[i].first,
                                              pred_vec[i].second.first,
                                              pred_vec[i].second.second);
            c = L_S.get_C(pred_vec[i].first);
            std::vector<std::tuple<uint64_t, word_t, std::pair<uint64_t, uint64_t>>> subj_vec;
            L_S.all_active_s_values_in_range_test<word_t>(c + interval.first, c + interval.second,
                                                          D_array, (word_t) A.next(current_D, pred_vec[i].first, BWD),
                                                          subj_vec);

            //PART3: Map the range of each subject to the range of objects
            for (uint64_t j = 0; j < subj_vec.size(); j++) {
                push_merge_interval(ist_container, subj_vec[j]);
                if (get<0>(subj_vec[j]) == const_2 && A.atFinal(get<1>(subj_vec[j]), BWD)) {
                    return true;
                }
            }
        }
        return false;

    }

    void push_merge_interval(Container& ist_container,
                             const std::tuple<uint64_t, word_t, std::pair<uint64_t, uint64_t>> &element){

        /**
             * get<0>: symbol
             * get<1>: NFA states (D)
             * get<2>: [lb, rb]
             */

        const auto lb = L_P.get_C(get<0>(element));
        const auto rb = L_P.get_C(get<0>(element) + 1) - 1;

        if(ist_container.empty()){
            ist_container.push(interval_state_type{bwt_interval(lb, rb),get<1>(element)});
            return;
        }
        //The previous interval is contiguous to element's interval and
        //both have equal NFA state
        auto& last = last_element(ist_container);
        if(last.interval.right()+1 == lb && last.current_D == get<1>(element)){
            last.interval.set_right(rb);
        }else{
            ist_container.push(interval_state_type{bwt_interval(lb, rb),get<1>(element)});
        }
    }


    void push_merge_interval_stats(Container& ist_container,
                             const std::tuple<uint64_t, word_t, std::pair<uint64_t, uint64_t>> &element,
                             uint64_t &merged, uint64_t &total, uint64_t &max_size){

        /**
             * get<0>: symbol
             * get<1>: NFA states (D)
             * get<2>: [lb, rb]
             */

        const auto lb = L_P.get_C(get<0>(element));
        const auto rb = L_P.get_C(get<0>(element) + 1) - 1;

        if(ist_container.empty()){
            ist_container.push(interval_state_type{bwt_interval(lb, rb),get<1>(element)});
            ++total;
            max_size = 1;
            return;
        }
        //The previous interval is contiguous to element's interval and
        //both have equal NFA state
        auto& last = last_element(ist_container);
        if(last.interval.right()+1 == lb && last.current_D == get<1>(element)){
            last.interval.set_right(rb);
            ++merged;
        }else{
            ist_container.push(interval_state_type{bwt_interval(lb, rb),get<1>(element)});
        }
        ++total;
        if(ist_container.size() > max_size){
            max_size = ist_container.size();
        }
    }









    bool rpq_var_to_var_obtain_o_with_bound(RpqAutomata &A, std::vector<uint64_t> &object_vector,
                                            std::vector<word_t> &B_array, high_resolution_clock::time_point start,
                                            uint64_t bound) {

        high_resolution_clock::time_point stop;
        double total_time = 0.0;
        duration<double> time_span;

        initializable_array<word_t> D_array(4 * (max_O + 1), 0);

        bool time_out = false;
        Container ist_container; //contains intervals with NFA states
        ist_container.push(interval_state_type{bwt_interval(1, nTriples), (word_t) A.getFinalStates()});
        while (!time_out && object_vector.size() < bound && !ist_container.empty()) {
            auto ist_top = first_element(ist_container);
            ist_container.pop();
            next_step_objects_with_bound(A, B_array, D_array, ist_top.current_D, ist_top.interval,
                                         ist_container, object_vector, bound);
            stop = high_resolution_clock::now();
            time_span = duration_cast<microseconds>(stop - start);
            total_time = time_span.count();
            if (total_time > TIME_OUT) time_out = true;
        }
    };

    bool rpq_var_to_var_obtain_o(RpqAutomata &A,
                                 std::vector<uint64_t> &object_vector,
                                 std::vector<word_t> &B_array,
                                 high_resolution_clock::time_point start) {

        high_resolution_clock::time_point stop;
        double total_time = 0.0;
        duration<double> time_span;

        initializable_array<word_t> D_array(4 * (max_O + 1), 0);

        bool time_out = false;
        Container ist_container; //contains intervals with NFA states
        ist_container.push(interval_state_type{bwt_interval(1, nTriples), (word_t) A.getFinalStates()});
        while (!time_out && !ist_container.empty()) {
            auto ist_top = first_element(ist_container);
            ist_container.pop();
            next_step_objects(A, B_array, D_array, ist_top.current_D,
                              ist_top.interval,ist_container, object_vector);
            stop = high_resolution_clock::now();
            time_span = duration_cast<microseconds>(stop - start);
            total_time = time_span.count();
            if (total_time > TIME_OUT) time_out = true;
        }
    };


    void _rpq_const_s_to_var_o_with_bound(RpqAutomata &A,
                               std::unordered_map<std::string, uint64_t> &predicates_map,
                               std::vector<word_t> &B_array,
                               uint64_t initial_object,
                               std::vector<std::pair<uint64_t, uint64_t>> &solutions,
                               bool const_to_var,
                               high_resolution_clock::time_point start,
                               uint64_t bound) {

        high_resolution_clock::time_point stop;
        double total_time = 0.0;
        duration<double> time_span;

        word_t current_D;  // palabra de maquina D, con los estados activos
        initializable_array<word_t> D_array(4 * (max_O + 1), 0);

        Container ist_container; //contains intervals with NFA states
        current_D = (word_t) A.getFinalStates();

        if (A.atFinal(current_D, BWD)) {
            solutions.emplace_back(initial_object, initial_object);
        }

        bool time_out = false;
        ist_container.push(interval_state_type{bwt_interval(L_P.get_C(initial_object),
                                                            L_P.get_C(initial_object + 1) - 1), current_D});
        while (!time_out && solutions.size() < bound && !ist_container.empty()) {
            auto ist_top = first_element(ist_container);
            ist_container.pop();
            next_step_const_to_var_with_bound(A, B_array, D_array, ist_top.current_D, ist_top.interval,
                                         initial_object, ist_container, solutions, bound, const_to_var);
            stop = high_resolution_clock::now();
            time_span = duration_cast<microseconds>(stop - start);
            total_time = time_span.count();
            if (total_time > TIME_OUT) time_out = true;
        }
    };

    bool _rpq_const_s_to_var_o(RpqAutomata &A,
                               std::unordered_map<std::string, uint64_t> &predicates_map,
                               std::vector<word_t> &B_array,
                               uint64_t initial_object,
                               std::vector<std::pair<uint64_t, uint64_t>> &solutions,
                               bool const_to_var,
                               high_resolution_clock::time_point start) {

        high_resolution_clock::time_point stop;
        double total_time = 0.0;
        duration<double> time_span;

        word_t current_D;  // palabra de maquina D, con los estados activos
        initializable_array<word_t> D_array(4 * (max_O + 1), 0);

        Container ist_container; //contains intervals with NFA states
        current_D = (word_t) A.getFinalStates();
        ist_container.push(interval_state_type{bwt_interval(L_P.get_C(initial_object),
                                                            L_P.get_C(initial_object + 1) - 1), current_D});
        if (A.atFinal(current_D, BWD)) {
            solutions.emplace_back(initial_object, initial_object);
        }

        while (!ist_container.empty()) {
            auto ist_top = first_element(ist_container);
            ist_container.pop();
            next_step_const_to_var_verbose(A, B_array, D_array, ist_top.current_D, ist_top.interval,
                                initial_object, ist_container, solutions, const_to_var);

            /*next_step_const_to_var(A, B_array, D_array, ist_top.current_D, ist_top.interval,
                                           initial_object, ist_container, solutions, const_to_var);*/
            stop = high_resolution_clock::now();
            time_span = duration_cast<microseconds>(stop - start);
            total_time = time_span.count();
            if (total_time > TIME_OUT) {
                //std::cout << "Time: " << total_time << std::endl;
                return true;
            }
        }
        return false;
    };

    bool _rpq_const_s_to_const_o(RpqAutomata &A,
                               std::unordered_map<std::string, uint64_t> &predicates_map,
                               std::vector<word_t> &B_array,
                               uint64_t const_1,
                               uint64_t const_2,
                               std::vector<std::pair<uint64_t, uint64_t>> &solutions,
                               bool const1_is_subj,
                               high_resolution_clock::time_point start) {

        high_resolution_clock::time_point stop;
        double total_time = 0.0;
        duration<double> time_span;

        word_t current_D;  // palabra de maquina D, con los estados activos
        initializable_array<word_t> D_array(4 * (max_O + 1), 0);

        Container ist_container; //contains intervals with NFA states
        current_D = (word_t) A.getFinalStates();
        ist_container.push(interval_state_type{bwt_interval(L_P.get_C(const_1),
                                                            L_P.get_C(const_1 + 1) - 1), current_D});
        if (const_1 == const_2 && A.atFinal(current_D, BWD)) {
            solutions.emplace_back(const_1, const_1);
            return false;
        }

        bool found = false;
        while (!ist_container.empty()) {
            auto ist_top = first_element(ist_container);
            ist_container.pop();
            found = next_step_const_to_const(A, B_array, D_array, ist_top.current_D, ist_top.interval,
                                        const_1, const_2, ist_container);
            stop = high_resolution_clock::now();
            time_span = duration_cast<microseconds>(stop - start);
            total_time = time_span.count();
            if (total_time > TIME_OUT) {
                //std::cout << "Time: " << total_time << std::endl;
                return true;
            }
            if(found){
                if(const1_is_subj){
                    solutions.emplace_back(const_1, const_2);
                }else{
                    solutions.emplace_back(const_2, const_1);
                }
                return false;
            }
        }
        return false;
    };

public:
    void or_query_var_to_var(const std::string &rpq, uint64_t n_or,
                             uint64_t bound, unordered_map<std::string, uint64_t> &predicates_map,
                             std::vector<std::pair<uint64_t, uint64_t>> &output) {
        uint64_t i, i1, pred, /*pred_1, pred_2,*/ k, start, a, max_pos;
        bool negated_pred/*, negated_p1, negated_p2*/;
        std::vector<uint64_t> values_x, values_y;
        std::pair<uint64_t, uint64_t> I_S;
        uint64_t c, object;
        std::unordered_set<std::pair<uint64_t, uint64_t>> o_set;
        std::pair<std::unordered_set<std::pair<uint64_t, uint64_t>>::iterator, bool> ret;
        std::vector<int64_t> pred_v;
        std::vector<bool> negated_pred_v;
        std::vector<std::pair<uint64_t, uint64_t>> I_S_v;

        for (i1 = k = 0; k < n_or; k++) {
            negated_pred_v.push_back(rpq.at(i1 + 1) == '%');
            start = i1;
            while (rpq.at(i1) != '>') ++i1;
            if (negated_pred_v[k])
                pred_v.push_back(predicates_map["<" + rpq.substr(start + 2, i1 - start - 1)]);
            else
                pred_v.push_back(predicates_map[rpq.substr(start, i1 - start + 1)]);

            I_S_v.push_back(std::pair<uint64_t, uint64_t>(
                    L_S.get_C(pred_v[k] + real_max_P),
                    L_S.get_C(pred_v[k] + real_max_P + 1) - 1
            ));
            i1 += 2;
        }

        for (i1 = k = 0; k < n_or and output.size() < bound; k++) {
            for (a = 0; a < n_or and pred_v[a] == -1; ++a);
            max_pos = a;
            for (++a; a < n_or; ++a) {
                if (pred_v[a] == -1) continue;
                if (I_S_v[a].second - I_S_v[a].first > I_S_v[max_pos].second - I_S_v[max_pos].first)
                    max_pos = a;
            }
            pred = pred_v[max_pos];
            I_S = I_S_v[max_pos];
            negated_pred = negated_pred_v[max_pos];
            pred_v[max_pos] = -1;
            // Now extract all ?x values
            L_S.all_values_in_range_bounded(I_S.first, I_S.second, values_x, bound - output.size());
            // For each ?x obtained, search for ?y p ?x using backward search
            for (i = 0; i < values_x.size() and output.size() < bound; i++) {
                object = values_x[i];
                I_S = L_P.backward_step(L_P.get_C(object), L_P.get_C(object + 1) - 1, pred);
                c = L_S.get_C(pred);
                I_S.first += c;
                I_S.second += c;
                values_y.clear();
                L_S.all_values_in_range_bounded(I_S.first, I_S.second, values_y, bound - output.size());
                if (negated_pred) {
                    for (uint64_t j = 0; j < values_y.size() and output.size() < bound; ++j) {
                        ret = o_set.insert(std::pair<uint64_t, uint64_t>(values_y[j], object));
                        if (ret.second == true)
                            output.push_back(std::pair<uint64_t, uint64_t>(values_y[j], object));
                    }
                } else {
                    for (uint64_t j = 0; j < values_y.size() and output.size() < bound; ++j) {
                        ret = o_set.insert(std::pair<uint64_t, uint64_t>(object, values_y[j]));
                        if (ret.second == true)
                            output.push_back(std::pair<uint64_t, uint64_t>(object, values_y[j]));
                    }
                }
            }

        }
    };

    void or_query_var_to_var(const std::string &rpq, uint64_t n_or,
                             unordered_map<std::string, uint64_t> &predicates_map,
                             std::vector<std::pair<uint64_t, uint64_t>> &output) {
        uint64_t i, i1, pred, /*pred_1, pred_2,*/ k, start, a, max_pos;
        bool negated_pred/*, negated_p1, negated_p2*/;
        std::vector<uint64_t> values_x, values_y;
        std::pair<uint64_t, uint64_t> I_S;
        uint64_t c, object;
        std::unordered_set<std::pair<uint64_t, uint64_t>> o_set;
        std::pair<std::unordered_set<std::pair<uint64_t, uint64_t>>::iterator, bool> ret;
        std::vector<int64_t> pred_v;
        std::vector<bool> negated_pred_v;
        std::vector<std::pair<uint64_t, uint64_t>> I_S_v;

        auto start_time = high_resolution_clock::now();
        for (i1 = k = 0; k < n_or; k++) {
            negated_pred_v.push_back(rpq.at(i1 + 1) == '%');
            start = i1;
            while (rpq.at(i1) != '>') ++i1;
            if (negated_pred_v[k])
                pred_v.push_back(predicates_map["<" + rpq.substr(start + 2, i1 - start - 1)]);
            else
                pred_v.push_back(predicates_map[rpq.substr(start, i1 - start + 1)]);

            I_S_v.push_back(std::pair<uint64_t, uint64_t>(
                    L_S.get_C(pred_v[k] + real_max_P),
                    L_S.get_C(pred_v[k] + real_max_P + 1) - 1
            ));
            i1 += 2;
        }


        auto stop = high_resolution_clock::now();
        auto total_time = duration_cast<seconds>(stop - start_time).count();
        bool time_out = (total_time > TIME_OUT);
        for (i1 = k = 0; !time_out && k < n_or; k++) {
            for (a = 0; a < n_or and pred_v[a] == -1; ++a);
            max_pos = a;
            for (++a; a < n_or; ++a) {
                if (pred_v[a] == -1) continue;
                if (I_S_v[a].second - I_S_v[a].first > I_S_v[max_pos].second - I_S_v[max_pos].first)
                    max_pos = a;
            }
            pred = pred_v[max_pos];
            I_S = I_S_v[max_pos];
            negated_pred = negated_pred_v[max_pos];
            pred_v[max_pos] = -1;
            // Now extract all ?x values
            values_x = L_S.all_values_in_range(I_S.first, I_S.second);
            stop = high_resolution_clock::now();
            total_time = duration_cast<seconds>(stop - start_time).count();
            time_out = (total_time > TIME_OUT);
            // For each ?x obtained, search for ?y p ?x using backward search
            for (i = 0; i < values_x.size() && !time_out; i++) {
                object = values_x[i];
                I_S = L_P.backward_step(L_P.get_C(object), L_P.get_C(object + 1) - 1, pred);
                c = L_S.get_C(pred);
                I_S.first += c;
                I_S.second += c;
                values_y.clear();
                values_y = L_S.all_values_in_range(I_S.first, I_S.second);
                if (negated_pred) {
                    for (uint64_t j = 0; !time_out && j < values_y.size(); ++j) {
                        ret = o_set.insert(std::pair<uint64_t, uint64_t>(values_y[j], object));
                        if (ret.second == true)
                            output.push_back(std::pair<uint64_t, uint64_t>(values_y[j], object));
                        stop = high_resolution_clock::now();
                        total_time = duration_cast<seconds>(stop - start_time).count();
                        time_out = (total_time > TIME_OUT);
                    }
                } else {
                    for (uint64_t j = 0; !time_out && j < values_y.size(); ++j) {
                        ret = o_set.insert(std::pair<uint64_t, uint64_t>(object, values_y[j]));
                        if (ret.second == true)
                            output.push_back(std::pair<uint64_t, uint64_t>(object, values_y[j]));
                        stop = high_resolution_clock::now();
                        total_time = duration_cast<seconds>(stop - start_time).count();
                        time_out = (total_time > TIME_OUT);
                    }
                }
            }

        }
    };


    // Solves single-predicate queries (with no operator, except ^)
    void single_predicate_query(uint64_t predicate, uint64_t object, uint64_t query_type,
                                bool is_negated, uint64_t bound,
                                std::vector<std::pair<uint64_t, uint64_t>> &output) {
        if (query_type == VAR_TO_CONST or query_type == CONST_TO_VAR) {
            std::pair<uint64_t, uint64_t> I_S = L_P.backward_step(L_P.get_C(object), L_P.get_C(object + 1) - 1,
                                                                  predicate);
            uint64_t c = L_S.get_C(predicate);
            I_S.first += c;
            I_S.second += c;

            std::vector<uint64_t> values_in_I_S;
            L_S.all_values_in_range_bounded(I_S.first, I_S.second, values_in_I_S, bound);

            if (VAR_TO_CONST and !is_negated)
                for (uint64_t i = 0; i < values_in_I_S.size(); ++i)
                    output.push_back(std::pair<uint64_t, uint64_t>(values_in_I_S[i], object));
            else
                for (uint64_t i = 0; i < values_in_I_S.size(); ++i)
                    output.push_back(std::pair<uint64_t, uint64_t>(object, values_in_I_S[i]));

        } else {
            // asumo que siempre envio el predicado original, sin negar, el var-to-var tiene que hacer eso
            std::pair<uint64_t, uint64_t> I_S = std::pair<uint64_t, uint64_t>(
                    L_S.get_C(predicate + real_max_P),
                    L_S.get_C(predicate + real_max_P + 1) - 1
            );

            // Now extract all ?x values
            std::vector<uint64_t> values_x, values_y;
            L_S.all_values_in_range_bounded(I_S.first, I_S.second, values_x, bound);
            uint64_t c, bound_aux = bound;
            // For each ?x obtained, search for ?y p ?x using backward search
            for (uint64_t i = 0; i < values_x.size() and output.size() < bound; i++) {
                object = values_x[i];
                I_S = L_P.backward_step(L_P.get_C(object), L_P.get_C(object + 1) - 1, predicate);
                c = L_S.get_C(predicate);
                I_S.first += c;
                I_S.second += c;
                values_y.clear();
                L_S.all_values_in_range_bounded(I_S.first, I_S.second, values_y, bound_aux);
                bound_aux -= values_y.size();
                for (uint64_t j = 0; j < values_y.size() and output.size() < bound; ++j)
                    output.push_back(std::pair<uint64_t, uint64_t>(object, values_y[j]));
            }
        }
    };

    // Solves single-predicate queries (with no operator, except ^)
    void single_predicate_query(uint64_t predicate, uint64_t object, uint64_t query_type,
                                bool is_negated,
                                std::vector<std::pair<uint64_t, uint64_t>> &output,
                                high_resolution_clock::time_point start) {
        if (query_type == VAR_TO_CONST or query_type == CONST_TO_VAR) {
            std::pair<uint64_t, uint64_t> I_S = L_P.backward_step(L_P.get_C(object), L_P.get_C(object + 1) - 1,
                                                                  predicate);
            uint64_t c = L_S.get_C(predicate);
            I_S.first += c;
            I_S.second += c;

            std::vector<uint64_t> values_in_I_S = L_S.all_values_in_range(I_S.first, I_S.second);

            if (VAR_TO_CONST and !is_negated) {
                auto stop = high_resolution_clock::now();
                auto total_time = duration_cast<seconds>(stop - start).count();
                bool time_out = (total_time > TIME_OUT);
                for (uint64_t i = 0; !time_out && i < values_in_I_S.size(); ++i) {
                    output.push_back(std::pair<uint64_t, uint64_t>(values_in_I_S[i], object));
                    stop = high_resolution_clock::now();
                    total_time = duration_cast<seconds>(stop - start).count();
                    time_out = (total_time > TIME_OUT);
                }
            }else{
                auto stop = high_resolution_clock::now();
                auto total_time = duration_cast<seconds>(stop - start).count();
                bool time_out = (total_time > TIME_OUT);
                for (uint64_t i = 0; !time_out && i < values_in_I_S.size(); ++i) {
                    output.push_back(std::pair<uint64_t, uint64_t>(object, values_in_I_S[i]));
                    stop = high_resolution_clock::now();
                    total_time = duration_cast<seconds>(stop - start).count();
                    time_out = (total_time > TIME_OUT);
                }

            }

        } else {
            // asumo que siempre envio el predicado original, sin negar, el var-to-var tiene que hacer eso
            std::pair<uint64_t, uint64_t> I_S = std::pair<uint64_t, uint64_t>(
                    L_S.get_C(predicate + real_max_P),
                    L_S.get_C(predicate + real_max_P + 1) - 1
            );

            // Now extract all ?x values
            std::vector<uint64_t> values_x, values_y;
            values_x = L_S.all_values_in_range(I_S.first, I_S.second);
            uint64_t c;
            auto stop = high_resolution_clock::now();
            auto total_time = duration_cast<seconds>(stop - start).count();
            bool time_out = (total_time > TIME_OUT);
            // For each ?x obtained, search for ?y p ?x using backward search
            for (uint64_t i = 0; !time_out && i < values_x.size(); i++) {
                object = values_x[i];
                I_S = L_P.backward_step(L_P.get_C(object), L_P.get_C(object + 1) - 1, predicate);
                c = L_S.get_C(predicate);
                I_S.first += c;
                I_S.second += c;
                values_y.clear();
                values_y = L_S.all_values_in_range(I_S.first, I_S.second);
                for (uint64_t j = 0; !time_out && j < values_y.size() ; ++j){
                    output.push_back(std::pair<uint64_t, uint64_t>(object, values_y[j]));
                    stop = high_resolution_clock::now();
                    total_time = duration_cast<seconds>(stop - start).count();
                    time_out = (total_time > TIME_OUT);
                }
            }
        }
    };

    void path_query(const std::string &rpq, uint64_t object, uint64_t query_type, uint64_t bound,
                    unordered_map<std::string, uint64_t> &predicates_map,
                    std::vector<std::pair<uint64_t, uint64_t>> &output) {
        if (query_type == VAR_TO_CONST || query_type == CONST_TO_VAR) {
            // primero voy a asumir que los predicados no son negados
            uint64_t i, pred_1, pred_2;
            for (i = 0; rpq.at(i) != '>'; ++i);

            if (query_type == VAR_TO_CONST) {
                pred_1 = predicates_map[rpq.substr(0, i + 1)];
                pred_2 = real_max_P + predicates_map[rpq.substr(i + 1, rpq.size() - 1)];  // actually, ^pred_2
            } else {
                pred_1 = real_max_P + predicates_map[rpq.substr(i + 1, rpq.size() - 1)];  // actually, ^pred_2
                pred_2 = predicates_map[rpq.substr(0, i + 1)];
            }

            // first, backward search for object ^pred_2, starting from the L_P interval of object
            std::pair<uint64_t, uint64_t> Is_p2 = L_P.backward_step(L_P.get_C(object), L_P.get_C(object + 1) - 1,
                                                                    pred_2);
            uint64_t c = L_S.get_C(pred_2);
            Is_p2.first += c;
            Is_p2.second += c;

            std::pair<uint64_t, uint64_t>
                    Is_p1 = std::pair<uint64_t, uint64_t>(L_S.get_C(pred_1),
                                                          L_S.get_C(pred_1 + 1) - 1);
            std::vector<std::pair<uint64_t, uint64_t>> z_values;
            std::vector<std::array<uint64_t, 2ul>> ranges;

            ranges.push_back({Is_p1.first, Is_p1.second});
            ranges.push_back({Is_p2.first, Is_p2.second});
            z_values = L_S.intersect(ranges);

            std::vector<uint64_t> values_s;
            uint64_t bound_aux = bound;
            if (query_type == VAR_TO_CONST)
                pred_1 += real_max_P; // now it becomes ^pred_1
            else
                pred_1 = pred_2 + real_max_P;
            std::unordered_set<std::pair<uint64_t, uint64_t>> o_set;
            std::pair<std::unordered_set<std::pair<uint64_t, uint64_t>>::iterator, bool> ret;
            for (i = 0; i < z_values.size() and output.size() < bound; ++i) {
                Is_p1 = L_P.backward_step(L_P.get_C(z_values[i].first), L_P.get_C(z_values[i].first + 1) - 1, pred_1);
                c = L_S.get_C(pred_1);
                Is_p1.first += c;
                Is_p1.second += c;
                values_s.clear();
                L_S.all_values_in_range_bounded(Is_p1.first, Is_p1.second, values_s, bound - output.size());
                bound_aux -= values_s.size();
                for (uint64_t j = 0; j < values_s.size() and output.size() < bound; ++j) {
                    ret = o_set.insert(std::pair<uint64_t, uint64_t>(values_s[j], z_values[i].first));
                    if (ret.second == true)
                        output.push_back(std::pair<uint64_t, uint64_t>(values_s[j], z_values[i].first));
                }
            }
        } else {
            if (query_type == VAR_TO_VAR) {
                // primero voy a asumir que los predicados no son negados
                uint64_t i;
                bool negated_p1, negated_p2;

                for (i = 0; rpq.at(i) != '>'; ++i);

                negated_p1 = rpq.at(1) == '%';
                negated_p2 = rpq.at(i + 2) == '%';

                uint64_t pred_1, pred_2;
                // first, negates pred_1
                if (negated_p1) {
                    pred_1 = real_max_P + predicates_map["<" + rpq.substr(2, i - 1)];
                    //cout << "<"+ rpq.substr(2, i-1) << endl;
                } else
                    pred_1 = predicates_map[rpq.substr(0, i + 1)];
                // then, negates pred_2
                if (negated_p2)
                    pred_2 = predicates_map["<" + rpq.substr(i + 3, rpq.size() - 1)];
                else
                    pred_2 = real_max_P + predicates_map[rpq.substr(i + 1, rpq.size() - 1)];  // actually, ^pred_2

                std::pair<uint64_t, uint64_t>
                        Is_p1 = std::pair<uint64_t, uint64_t>(L_S.get_C(pred_1),
                                                              L_S.get_C(pred_1 + 1) - 1),
                        Is_p2 = std::pair<uint64_t, uint64_t>(L_S.get_C(pred_2),
                                                              L_S.get_C(pred_2 + 1) - 1);

                std::vector<std::pair<uint64_t, uint64_t>> z_values;
                std::vector<std::array<uint64_t, 2ul>> ranges;

                ranges.push_back({Is_p1.first, Is_p1.second});
                ranges.push_back({Is_p2.first, Is_p2.second});
                z_values = L_S.intersect(ranges);

                std::vector<uint64_t> values_s1, values_s2;
                uint64_t c, z;

                pred_1 = pred_1 + ((negated_p1) ? -real_max_P : real_max_P);
                pred_2 = pred_2 + ((negated_p2) ? real_max_P : -real_max_P); // now it becomes pred_2

                std::unordered_set<std::pair<uint64_t, uint64_t>> o_set;
                std::pair<std::unordered_set<std::pair<uint64_t, uint64_t>>::iterator, bool> ret;
                for (i = 0; i < z_values.size() and output.size() < bound; ++i) {
                    z = z_values[i].first;
                    Is_p1 = L_P.backward_step(L_P.get_C(z), L_P.get_C(z + 1) - 1, pred_1);
                    c = L_S.get_C(pred_1);
                    Is_p1.first += c;
                    Is_p1.second += c;
                    Is_p2 = L_P.backward_step(L_P.get_C(z), L_P.get_C(z + 1) - 1, pred_2);
                    c = L_S.get_C(pred_2);
                    Is_p2.first += c;
                    Is_p2.second += c;
                    values_s1.clear();
                    values_s2.clear();
                    L_S.all_values_in_range_bounded(Is_p1.first, Is_p1.second, values_s1, bound - output.size());
                    L_S.all_values_in_range_bounded(Is_p2.first, Is_p2.second, values_s2, bound - output.size());

                    for (uint64_t j = 0; j < values_s1.size() and output.size() < bound; ++j) {
                        for (uint64_t w = 0; w < values_s2.size() and output.size() < bound; ++w) {
                            ret = o_set.insert(std::pair<uint64_t, uint64_t>(values_s1[j], values_s2[w]));
                            if (ret.second == true)
                                output.push_back(std::pair<uint64_t, uint64_t>(values_s1[j], values_s2[w]));
                        }
                    }
                }
            }
        }
    };

    void path_query(const std::string &rpq, uint64_t object, uint64_t query_type,
                    unordered_map<std::string, uint64_t> &predicates_map,
                    std::vector<std::pair<uint64_t, uint64_t>> &output,
                    high_resolution_clock::time_point start) {

        if (query_type == VAR_TO_CONST || query_type == CONST_TO_VAR) {
            // primero voy a asumir que los predicados no son negados
            uint64_t i, pred_1, pred_2;
            for (i = 0; rpq.at(i) != '>'; ++i);

            if (query_type == VAR_TO_CONST) {
                pred_1 = predicates_map[rpq.substr(0, i + 1)];
                pred_2 = real_max_P + predicates_map[rpq.substr(i + 1, rpq.size() - 1)];  // actually, ^pred_2
            } else {
                pred_1 = real_max_P + predicates_map[rpq.substr(i + 1, rpq.size() - 1)];  // actually, ^pred_2
                pred_2 = predicates_map[rpq.substr(0, i + 1)];
            }

            // first, backward search for object ^pred_2, starting from the L_P interval of object
            std::pair<uint64_t, uint64_t> Is_p2 = L_P.backward_step(L_P.get_C(object), L_P.get_C(object + 1) - 1,
                                                                    pred_2);
            uint64_t c = L_S.get_C(pred_2);
            Is_p2.first += c;
            Is_p2.second += c;

            std::pair<uint64_t, uint64_t>
                    Is_p1 = std::pair<uint64_t, uint64_t>(L_S.get_C(pred_1),
                                                          L_S.get_C(pred_1 + 1) - 1);
            std::vector<std::pair<uint64_t, uint64_t>> z_values;
            std::vector<std::array<uint64_t, 2ul>> ranges;

            ranges.push_back({Is_p1.first, Is_p1.second});
            ranges.push_back({Is_p2.first, Is_p2.second});
            z_values = L_S.intersect(ranges);

            std::vector<uint64_t> values_s;
            if (query_type == VAR_TO_CONST)
                pred_1 += real_max_P; // now it becomes ^pred_1
            else
                pred_1 = pred_2 + real_max_P;
            std::unordered_set<std::pair<uint64_t, uint64_t>> o_set;
            std::pair<std::unordered_set<std::pair<uint64_t, uint64_t>>::iterator, bool> ret;

            auto stop = high_resolution_clock::now();
            auto total_time = duration_cast<seconds>(stop - start).count();
            auto time_out = (total_time > TIME_OUT);
            for (i = 0; !time_out && i < z_values.size(); ++i) {
                Is_p1 = L_P.backward_step(L_P.get_C(z_values[i].first), L_P.get_C(z_values[i].first + 1) - 1, pred_1);
                c = L_S.get_C(pred_1);
                Is_p1.first += c;
                Is_p1.second += c;
                values_s.clear();
                values_s = L_S.all_values_in_range(Is_p1.first, Is_p1.second);
                for (uint64_t j = 0; !time_out && j < values_s.size(); ++j) {
                    ret = o_set.insert(std::pair<uint64_t, uint64_t>(values_s[j], z_values[i].first));
                    if (ret.second == true)
                        output.push_back(std::pair<uint64_t, uint64_t>(values_s[j], z_values[i].first));
                    stop = high_resolution_clock::now();
                    total_time = duration_cast<seconds>(stop - start).count();
                    time_out = (total_time > TIME_OUT);
                }
            }
        } else {
            if (query_type == VAR_TO_VAR) {
                // primero voy a asumir que los predicados no son negados
                uint64_t i;
                bool negated_p1, negated_p2;

                for (i = 0; rpq.at(i) != '>'; ++i);

                negated_p1 = rpq.at(1) == '%';
                negated_p2 = rpq.at(i + 2) == '%';

                uint64_t pred_1, pred_2;
                // first, negates pred_1
                if (negated_p1) {
                    pred_1 = real_max_P + predicates_map["<" + rpq.substr(2, i - 1)];
                    //cout << "<"+ rpq.substr(2, i-1) << endl;
                } else
                    pred_1 = predicates_map[rpq.substr(0, i + 1)];
                // then, negates pred_2
                if (negated_p2)
                    pred_2 = predicates_map["<" + rpq.substr(i + 3, rpq.size() - 1)];
                else
                    pred_2 = real_max_P + predicates_map[rpq.substr(i + 1, rpq.size() - 1)];  // actually, ^pred_2

                std::pair<uint64_t, uint64_t>
                        Is_p1 = std::pair<uint64_t, uint64_t>(L_S.get_C(pred_1),
                                                              L_S.get_C(pred_1 + 1) - 1),
                        Is_p2 = std::pair<uint64_t, uint64_t>(L_S.get_C(pred_2),
                                                              L_S.get_C(pred_2 + 1) - 1);

                std::vector<std::pair<uint64_t, uint64_t>> z_values;
                std::vector<std::array<uint64_t, 2ul>> ranges;

                ranges.push_back({Is_p1.first, Is_p1.second});
                ranges.push_back({Is_p2.first, Is_p2.second});
                z_values = L_S.intersect(ranges);

                std::vector<uint64_t> values_s1, values_s2;
                uint64_t c, z;

                pred_1 = pred_1 + ((negated_p1) ? -real_max_P : real_max_P);
                pred_2 = pred_2 + ((negated_p2) ? real_max_P : -real_max_P); // now it becomes pred_2

                std::unordered_set<std::pair<uint64_t, uint64_t>> o_set;
                std::pair<std::unordered_set<std::pair<uint64_t, uint64_t>>::iterator, bool> ret;
                auto stop = high_resolution_clock::now();
                auto total_time = duration_cast<seconds>(stop - start).count();
                auto time_out = (total_time > TIME_OUT);
                for (i = 0; !time_out && i < z_values.size(); ++i) {
                    z = z_values[i].first;
                    Is_p1 = L_P.backward_step(L_P.get_C(z), L_P.get_C(z + 1) - 1, pred_1);
                    c = L_S.get_C(pred_1);
                    Is_p1.first += c;
                    Is_p1.second += c;
                    Is_p2 = L_P.backward_step(L_P.get_C(z), L_P.get_C(z + 1) - 1, pred_2);
                    c = L_S.get_C(pred_2);
                    Is_p2.first += c;
                    Is_p2.second += c;
                    values_s1.clear();
                    values_s2.clear();
                    values_s1 = L_S.all_values_in_range(Is_p1.first, Is_p1.second);
                    values_s2 = L_S.all_values_in_range(Is_p2.first, Is_p2.second);

                    for (uint64_t j = 0; !time_out && j < values_s1.size(); ++j) {
                        for (uint64_t w = 0; !time_out && w < values_s2.size(); ++w) {
                            ret = o_set.insert(std::pair<uint64_t, uint64_t>(values_s1[j], values_s2[w]));
                            if (ret.second == true)
                                output.push_back(std::pair<uint64_t, uint64_t>(values_s1[j], values_s2[w]));
                            stop = high_resolution_clock::now();
                            total_time = duration_cast<seconds>(stop - start).count();
                            time_out = (total_time > TIME_OUT);
                        }
                    }
                }
            }
        }
    };



    void rpq_const_s_to_var_o(const std::string &rpq,
                              unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                              std::vector<word_t> &B_array,
                              uint64_t initial_object,
                              std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                              uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path) {

        std::string query, str_aux;


        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = real_max_P + predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.
            high_resolution_clock::time_point start = high_resolution_clock::now();
            single_predicate_query(predicate, initial_object, CONST_TO_VAR, is_negated_pred,
                                   output_subjects, start);
            return;
        } else {
            if (is_a_path and n_operators == 1) {
                high_resolution_clock::time_point start = high_resolution_clock::now();
                path_query(rpq, initial_object, CONST_TO_VAR, predicates_map, output_subjects, start);
                return;
            }
        }

        int64_t iii = rpq.size() - 1;
        query = parse_reverse(rpq, iii, predicates_map, real_max_P);

        RpqAutomata A(query, predicates_map);

        std::unordered_map<uint64_t, uint64_t> m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++) {
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);
        }

        high_resolution_clock::time_point start;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        _rpq_const_s_to_var_o(A, predicates_map, B_array, initial_object, output_subjects, true, start);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++) {
            L_P.unmark<word_t>(it->first, B_array);
        }
    };


    void rpq_const_s_to_var_o(const std::string &rpq,
                              unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                              std::vector<word_t> &B_array,
                              uint64_t initial_object,
                              std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                              uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path,
                              uint64_t bound) {
        std::string query, str_aux;

        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = real_max_P + predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.

            single_predicate_query(predicate, initial_object, CONST_TO_VAR, is_negated_pred, bound,
                                   output_subjects);
            return;
        } else {
            if (is_a_path and n_operators == 1) {
                path_query(rpq, initial_object, CONST_TO_VAR, bound, predicates_map, output_subjects);
                return;
            }
        }

        int64_t iii = rpq.size() - 1;
        query = parse_reverse(rpq, iii, predicates_map, real_max_P);

        RpqAutomata A(query, predicates_map);

        std::unordered_map<uint64_t, uint64_t> m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++) {
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);
        }

        high_resolution_clock::time_point start;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        _rpq_const_s_to_var_o_with_bound(A, predicates_map, B_array, initial_object, output_subjects, true, start, bound);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++) {
            L_P.unmark<word_t>(it->first, B_array);
        }
    };

    void rpq_var_s_to_const_o(const std::string &rpq,
                              unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                              std::vector<word_t> &B_array,
                              uint64_t initial_object,
                              std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                              uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path) {
        std::string query, str_aux;

        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = real_max_P + predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.

            high_resolution_clock::time_point start = high_resolution_clock::now();
            single_predicate_query(predicate, initial_object, VAR_TO_CONST, is_negated_pred,
                                   output_subjects, start);
            return;

        } else {
            if (is_a_path and n_operators == 1) {
                high_resolution_clock::time_point start = high_resolution_clock::now();
                path_query(rpq, initial_object, VAR_TO_CONST, predicates_map, output_subjects, start);
                return;
            }
        }


        int64_t iii = 0;
        query = parse(rpq, iii, predicates_map, real_max_P);

        RpqAutomata A(query, predicates_map);

        // ToDo: actualmente asume que el arreglo B_array tiene espacio tambien para las hojas del WT
        // Se puede cambiar y reducir el espacio del arreglo a la mitad, manejando las hojas con
        // el unordered map
        std::unordered_map<uint64_t, uint64_t> m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++) {
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);
        }

        high_resolution_clock::time_point start;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        _rpq_const_s_to_var_o(A, predicates_map, B_array, initial_object, output_subjects, false, start);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++) {
            L_P.unmark<word_t>(it->first, B_array);
        }
    };


    void rpq_var_s_to_const_o(const std::string &rpq,
                              unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                              std::vector<word_t> &B_array,
                              uint64_t initial_object,
                              std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                              uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path,
                              uint64_t bound) {
        std::string query, str_aux;

        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = real_max_P + predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.

            single_predicate_query(predicate, initial_object, VAR_TO_CONST, is_negated_pred, bound,
                                   output_subjects);
            return;

        } else {
            if (is_a_path and n_operators == 1) {
                path_query(rpq, initial_object, VAR_TO_CONST, bound, predicates_map, output_subjects);
                return;
            }
        }


        int64_t iii = 0;
        query = parse(rpq, iii, predicates_map, real_max_P);

        RpqAutomata A(query, predicates_map);

        // ToDo: actualmente asume que el arreglo B_array tiene espacio tambien para las hojas del WT
        // Se puede cambiar y reducir el espacio del arreglo a la mitad, manejando las hojas con
        // el unordered map
        std::unordered_map<uint64_t, uint64_t> m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++) {
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);
        }

        high_resolution_clock::time_point start;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        _rpq_const_s_to_var_o_with_bound(A, predicates_map, B_array, initial_object, output_subjects, false, start, bound);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++) {
            L_P.unmark<word_t>(it->first, B_array);
        }
    };



    void rpq_var_s_to_var_o(const std::string &rpq,
                            unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                            std::vector<word_t> &B_array,
                            std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                            uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path) {
        std::string query, str_aux;

        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.
            high_resolution_clock::time_point start = high_resolution_clock::now();
            single_predicate_query(predicate, 0, VAR_TO_VAR, is_negated_pred,
                                   output_subjects, start);
            return;
        } else {
            if (is_a_path and n_operators == 1) {
                high_resolution_clock::time_point start = high_resolution_clock::now();
                path_query(rpq, 0, VAR_TO_VAR, predicates_map, output_subjects, start);
                return;
            }
        }

        // First o -> s

        int64_t i = 0;
        query = parse(rpq, i, predicates_map, real_max_P);
        RpqAutomata A(query, predicates_map);

        std::unordered_map<uint64_t, uint64_t> m = A.getB();
        high_resolution_clock::time_point start, stop;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        str_aux.clear();

        // then s -> o
        std::string query2;

        i = rpq.size() - 1;
        query2 = parse_reverse(rpq, i, predicates_map, real_max_P);

        RpqAutomata A2(query2, predicates_map);
        m = A2.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        std::vector<uint64_t> object_vector_2;

        rpq_var_to_var_obtain_o(A2, object_vector_2, B_array, start);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);

        m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        bool time_out = false;
        for (uint64_t i = 0; !time_out && i < object_vector_2.size(); i++)
            time_out = _rpq_const_s_to_var_o(A, predicates_map, B_array, object_vector_2[i], output_subjects, true, start);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);
    };

    void rpq_var_s_to_var_o(const std::string &rpq,
                            unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                            std::vector<word_t> &B_array,
                            std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                            uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path,
                            uint64_t bound) {
        std::string query, str_aux;

        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.
            single_predicate_query(predicate, 0, VAR_TO_VAR, is_negated_pred, bound,
                                   output_subjects);
            return;
        } else {
            if (is_a_path and n_operators == 1) {
                path_query(rpq, 0, VAR_TO_VAR, bound, predicates_map, output_subjects);
                return;
            }
        }

        // First o -> s

        int64_t i = 0;
        query = parse(rpq, i, predicates_map, real_max_P);
        RpqAutomata A(query, predicates_map);

        std::unordered_map<uint64_t, uint64_t> m = A.getB();
        high_resolution_clock::time_point start, stop;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        str_aux.clear();

        // then s -> o
        std::string query2;

        i = rpq.size() - 1;
        query2 = parse_reverse(rpq, i, predicates_map, real_max_P);

        RpqAutomata A2(query2, predicates_map);
        m = A2.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        std::vector<uint64_t> object_vector_2;

        rpq_var_to_var_obtain_o_with_bound(A2, object_vector_2, B_array, start, bound);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);

        m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        for (uint64_t i = 0; i < object_vector_2.size() and output_subjects.size() < bound; i++)
            _rpq_const_s_to_var_o_with_bound(A, predicates_map, B_array, object_vector_2[i], output_subjects, true, start, bound);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);
    };


    void rpq_var_to_var_min(const std::string &rpq,
                            unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                            std::vector<word_t> &B_array,
                            std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                            uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path,
                            uint64_t bound) {
        std::string query, str_aux;

        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.
            single_predicate_query(predicate, 0, VAR_TO_VAR, is_negated_pred, bound,
                                   output_subjects);
            return;
        } else {
            if (is_a_path and n_operators == 1) {
                path_query(rpq, 0, VAR_TO_VAR, bound, predicates_map, output_subjects);
                return;
            }
        }

        // First o -> s
        int64_t i = 0;
        query = parse(rpq, i, predicates_map, real_max_P);
        RpqAutomata A(query, predicates_map);

        // ToDo: actualmente asume que el arreglo B_array tiene espacio tambien para las hojas del WT
        // Se puede cambiar y reducir el espacio del arreglo a la mitad, manejando las hojas con
        // el unordered map
        std::unordered_map<uint64_t, uint64_t> m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        std::vector<uint64_t> object_vector_1;

        high_resolution_clock::time_point start, stop;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        rpq_var_to_var_obtain_o_with_bound(A, object_vector_1, B_array, start, bound);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);

        str_aux.clear();

        // then s -> o
        std::string query2;
        i = rpq.size() - 1;
        query2 = parse_reverse(rpq, i, predicates_map, real_max_P);
        RpqAutomata A2(query2, predicates_map);
        m = A2.getB();
        std::vector<uint64_t> object_vector_2;

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        rpq_var_to_var_obtain_o_with_bound(A2, object_vector_2, B_array, start, bound);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);

        if (object_vector_1.size() < object_vector_2.size()) {
            m = A2.getB();
            for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
                L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

            for (uint64_t i = 0; i < object_vector_1.size() and output_subjects.size() < bound; i++)
                _rpq_const_s_to_var_o_with_bound(A2, predicates_map, B_array, object_vector_1[i], output_subjects, true, start,
                                      bound);
        } else {
            m = A.getB();
            for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
                L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

            for (uint64_t i = 0; i < object_vector_2.size() and output_subjects.size() < bound; i++)
                _rpq_const_s_to_var_o_with_bound(A, predicates_map, B_array, object_vector_2[i], output_subjects, true, start,
                                      bound);
        }

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);
    };

    void rpq_var_to_var_min(const std::string &rpq,
                            unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                            std::vector<word_t> &B_array,
                            std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                            uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path) {
        std::string query, str_aux;

        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.

            high_resolution_clock::time_point start = high_resolution_clock::now();
            single_predicate_query(predicate, 0, VAR_TO_VAR, is_negated_pred,
                                   output_subjects, start);
            return;
        } else {
            if (is_a_path and n_operators == 1) {
                high_resolution_clock::time_point start = high_resolution_clock::now();
                path_query(rpq, 0, VAR_TO_VAR, predicates_map, output_subjects, start);
                return;
            }
        }

        // First o -> s
        int64_t i = 0;
        query = parse(rpq, i, predicates_map, real_max_P);
        RpqAutomata A(query, predicates_map);

        // ToDo: actualmente asume que el arreglo B_array tiene espacio tambien para las hojas del WT
        // Se puede cambiar y reducir el espacio del arreglo a la mitad, manejando las hojas con
        // el unordered map
        std::unordered_map<uint64_t, uint64_t> m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        std::vector<uint64_t> object_vector_1;

        high_resolution_clock::time_point start, stop;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        rpq_var_to_var_obtain_o(A, object_vector_1, B_array, start);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);

        str_aux.clear();

        // then s -> o
        std::string query2;
        i = rpq.size() - 1;
        query2 = parse_reverse(rpq, i, predicates_map, real_max_P);
        RpqAutomata A2(query2, predicates_map);
        m = A2.getB();
        std::vector<uint64_t> object_vector_2;

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        rpq_var_to_var_obtain_o(A2, object_vector_2, B_array, start);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);

        if (object_vector_1.size() < object_vector_2.size()) {
            m = A2.getB();
            for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
                L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

            bool time_out = false;
            for (uint64_t i = 0; !time_out && i < object_vector_1.size(); i++)
                time_out = _rpq_const_s_to_var_o(A2, predicates_map, B_array, object_vector_1[i], output_subjects, true, start);
        } else {
            m = A.getB();
            for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
                L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

            bool time_out = false;
            for (uint64_t i = 0; !time_out && i < object_vector_2.size(); i++)
                time_out = _rpq_const_s_to_var_o(A, predicates_map, B_array, object_vector_2[i], output_subjects, true, start);
        }

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);
    };


    void rpq_var_to_var_os(const std::string &rpq,
                           unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                           std::vector<word_t> &B_array,
                           std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                           uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path,
                           uint64_t bound) {
        std::string query, str_aux;

        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.
            single_predicate_query(predicate, 0, VAR_TO_VAR, is_negated_pred, bound,
                                   output_subjects);
            return;
        } else {
            if (is_a_path and n_operators == 1) {
                path_query(rpq, 0, VAR_TO_VAR, bound, predicates_map, output_subjects);
                return;
            }
        }

        int64_t i = 0;
        query = parse(rpq, i, predicates_map, real_max_P);
        RpqAutomata A(query, predicates_map);

        // ToDo: actualmente asume que el arreglo B_array tiene espacio tambien para las hojas del WT
        // Se puede cambiar y reducir el espacio del arreglo a la mitad, manejando las hojas con
        // el unordered map
        std::unordered_map<uint64_t, uint64_t> m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        std::vector<uint64_t> object_vector_1;

        high_resolution_clock::time_point start, stop;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        rpq_var_to_var_obtain_o_with_bound(A, object_vector_1, B_array, start, bound);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);

        str_aux.clear();

        // then s -> o
        std::string query2;

        i = rpq.size() - 1;
        query2 = parse_reverse(rpq, i, predicates_map, real_max_P);

        RpqAutomata A2(query2, predicates_map);

        m = A2.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        for (uint64_t i = 0; i < object_vector_1.size() and output_subjects.size() < bound; i++)
            _rpq_const_s_to_var_o_with_bound(A2, predicates_map, B_array, object_vector_1[i], output_subjects, true, start, bound);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);
    };

    void rpq_var_to_var_os(const std::string &rpq,
                           unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                           std::vector<word_t> &B_array,
                           std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                           uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path) {
        std::string query, str_aux;

        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.
            high_resolution_clock::time_point start = high_resolution_clock::now();
            single_predicate_query(predicate, 0, VAR_TO_VAR, is_negated_pred,
                                   output_subjects, start);
            return;
        } else {
            if (is_a_path and n_operators == 1) {
                high_resolution_clock::time_point start = high_resolution_clock::now();
                path_query(rpq, 0, VAR_TO_VAR, predicates_map, output_subjects, start);
                return;
            }
        }

        int64_t i = 0;
        query = parse(rpq, i, predicates_map, real_max_P);
        RpqAutomata A(query, predicates_map);

        // ToDo: actualmente asume que el arreglo B_array tiene espacio tambien para las hojas del WT
        // Se puede cambiar y reducir el espacio del arreglo a la mitad, manejando las hojas con
        // el unordered map
        std::unordered_map<uint64_t, uint64_t> m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        std::vector<uint64_t> object_vector_1;

        high_resolution_clock::time_point start, stop;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        rpq_var_to_var_obtain_o(A, object_vector_1, B_array, start);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);


        str_aux.clear();

        // then s -> o
        std::string query2;

        i = rpq.size() - 1;
        query2 = parse_reverse(rpq, i, predicates_map, real_max_P);

        RpqAutomata A2(query2, predicates_map);

        m = A2.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        bool time_out = false;
        for (uint64_t i = 0; !time_out && i < object_vector_1.size(); i++)
            time_out = _rpq_const_s_to_var_o(A2, predicates_map, B_array, object_vector_1[i], output_subjects, true, start);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);
    };


    void rpq_var_to_var_so(const std::string &rpq,
                           unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                           std::vector<word_t> &B_array,
                           std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                           uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path,
                           uint64_t bound) {
        std::string query, str_aux;

        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.
            single_predicate_query(predicate, 0, VAR_TO_VAR, is_negated_pred, bound,
                                   output_subjects);
            return;
        } else {
            if (is_a_path and n_operators == 1) {
                path_query(rpq, 0, VAR_TO_VAR, bound, predicates_map, output_subjects);
                return;
            }
        }

        int64_t i = 0;
        query = parse(rpq, i, predicates_map, real_max_P);
        RpqAutomata A(query, predicates_map);

        // ToDo: actualmente asume que el arreglo B_array tiene espacio tambien para las hojas del WT
        // Se puede cambiar y reducir el espacio del arreglo a la mitad, manejando las hojas con
        // el unordered map
        std::unordered_map<uint64_t, uint64_t> m;
        high_resolution_clock::time_point start, stop;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        str_aux.clear();

        std::string query2;

        i = rpq.size() - 1;
        query2 = parse_reverse(rpq, i, predicates_map, real_max_P);

        RpqAutomata A2(query2, predicates_map);
        m = A2.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        std::vector<uint64_t> object_vector_2;

        rpq_var_to_var_obtain_o_with_bound(A2, object_vector_2, B_array, start, bound);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);

        m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        for (uint64_t i = 0; i < object_vector_2.size() and output_subjects.size() < bound; i++)
            _rpq_const_s_to_var_o_with_bound(A, predicates_map, B_array, object_vector_2[i], output_subjects, true, start, bound);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);

    };

    void rpq_var_to_var_so(const std::string &rpq,
                           unordered_map<std::string, uint64_t> &predicates_map,  // ToDo: esto debería ser una variable miembro de la clase
                           std::vector<word_t> &B_array,
                           std::vector<std::pair<uint64_t, uint64_t>> &output_subjects,
                           uint64_t n_predicates, bool is_negated_pred, uint64_t n_operators, bool is_a_path) {
        std::string query, str_aux;

        if (n_predicates == 1 and n_operators == 0) {
            uint64_t predicate;
            if (is_negated_pred)
                predicate = predicates_map["<" + rpq.substr(2, rpq.size() - 1)];
            else
                predicate = predicates_map[rpq];

            // cuidado, lo anterior asume que los negados han sido manipulados desde afuera, lo cual es cierto en la manera que los estoy escribiendo en el log que manejo, pero hay que hacerlo de una forma mas general.
            high_resolution_clock::time_point start = high_resolution_clock::now();
            single_predicate_query(predicate, 0, VAR_TO_VAR, is_negated_pred,
                                   output_subjects, start);
            return;
        } else {
            if (is_a_path and n_operators == 1) {
                high_resolution_clock::time_point start = high_resolution_clock::now();
                path_query(rpq, 0, VAR_TO_VAR, predicates_map, output_subjects, start);
                return;
            }
        }

        int64_t i = 0;
        query = parse(rpq, i, predicates_map, real_max_P);
        RpqAutomata A(query, predicates_map);

        // ToDo: actualmente asume que el arreglo B_array tiene espacio tambien para las hojas del WT
        // Se puede cambiar y reducir el espacio del arreglo a la mitad, manejando las hojas con
        // el unordered map
        std::unordered_map<uint64_t, uint64_t> m;
        high_resolution_clock::time_point start, stop;
        double total_time = 0.0;
        duration<double> time_span;
        start = high_resolution_clock::now();

        str_aux.clear();

        std::string query2;

        i = rpq.size() - 1;
        query2 = parse_reverse(rpq, i, predicates_map, real_max_P);

        RpqAutomata A2(query2, predicates_map);
        m = A2.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);

        std::vector<uint64_t> object_vector_2;

        rpq_var_to_var_obtain_o(A2, object_vector_2, B_array, start);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);



        m = A.getB();
        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.mark<word_t>(it->first, B_array, (word_t) it->second);


        bool time_out = false;
        for (uint64_t i = 0; !time_out && i < object_vector_2.size(); i++)
            time_out = _rpq_const_s_to_var_o(A, predicates_map, B_array, object_vector_2[i], output_subjects, true, start);

        for (std::unordered_map<uint64_t, uint64_t>::iterator it = m.begin(); it != m.end(); it++)
            L_P.unmark<word_t>(it->first, B_array);

    };


};


typedef ring_rpq<std::stack<interval_state_type>> ring_rpq_basic_dfs;
typedef ring_rpq<std::queue<interval_state_type>> ring_rpq_basic_bfs;

#endif
