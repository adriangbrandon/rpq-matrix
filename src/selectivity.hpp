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
// Created by Adrián on 8/6/22.
//

#ifndef RING_RPQ_SELECTIVITY_HPP
#define RING_RPQ_SELECTIVITY_HPP

#include "bwt-C-nose.hpp"
#include "bwt_interval.hpp"
#include "RpqTree.hpp"


namespace selectivity {

    enum split_type {none, target, source, intersect};
    enum query_type {var_var, const_var, var_const};

    inline uint64_t reverse(uint64_t id, uint64_t maxP){
        return (id > maxP) ? id - maxP : id + maxP;
    }

    inline uint64_t distinct_values(const uint64_t l, const uint64_t r, bwt_type &wt_pred_s){
        return  wt_pred_s.count_range_search_2d(l, r, 0, l-1);
    }


    /*inline uint64_t distinct_values(const uint64_t l, const uint64_t r, bwt_nose &L_S){
        return  L_S.count_distinct_values(l, r);
    }*/

    struct info {
        double weight;
        split_type  split;
        bool first_left;
    };

    struct info_preds {
        double weight;
        split_type  split;
        uint64_t mand_pred_left;
        uint64_t mand_pred_right;
    };

    class h_distinct {
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<double> m_l;
        std::vector<double> m_r;
        uint64_t m_max_p;
        uint64_t m_sigma;

    public:



        h_distinct(const std::vector<PairPredPos> &preds,
                   bwt_nose &L_S, bwt_type &wt_pred_s,
                        uint64_t maxP, uint64_t sigma) {

            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            for (const auto &pair : preds) {

                auto e_d = L_S.get_C(pair.pred_e + 1) - 1;
                auto b_d = L_S.get_C(pair.pred_e);
                auto v_target = distinct_values(b_d, e_d, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.pred_b, m_max_p);
                auto e_r = L_S.get_C(rev_id + 1) - 1;
                auto b_r = L_S.get_C(rev_id);
                auto v_source = distinct_values(b_r, e_r, wt_pred_s);
                m_s.push_back(v_source);
            }
            auto s = m_s.size();
            m_r.resize(s);
            //m_r[s-1]=1;
            m_l.resize(s);
            // m_l[0]=1;
            for(uint64_t i = 0; i < s; ++i){
                m_r[i] = m_t[i] / (double) m_s[i];
                m_l[i] = m_s[i] / (double) m_t[i];
            }
        }

        info simple(const uint64_t ith){
            info res;
            //double a;
            if(m_s[ith] < m_t[ith]){
                res.split = source;
                res.weight = m_s[ith] / (double) (m_sigma);
                double w_r = 1, w_l = 1;
                for(uint64_t i = ith; i < m_r.size(); ++i){
                    w_r *= m_r[i];
                }
                for(int64_t i = ith-1; i >= 0; --i){
                    w_l *= m_l[i];
                }
                res.first_left = (w_r >= w_l);
            }else{
                res.split = target;
                res.weight = m_t[ith] / (double) (m_sigma);
                double w_r = 1, w_l = 1;
                for(uint64_t i = ith+1; i < m_r.size(); ++i){
                    w_r *= m_r[i];
                }
                for(int64_t i = ith; i >= 0; --i){
                    w_l *= m_l[i];
                }
                res.first_left = (w_r >= w_l);
            }


            return res;
        }

        info intersection(const uint64_t ith){
            info res;
            res.split = intersect;
            res.weight = (double) (m_s[ith+1] * m_t[ith]) / (double) (m_sigma * m_sigma);
            double w_r = 1, w_l = 1;
            for(uint64_t i = ith+1; i < m_r.size(); ++i){
                w_r *= m_r[i];
            }
            for(int64_t i = ith; i >= 0; --i){
                w_l *= m_l[i];
            }
            res.first_left = (w_r >= w_l);
            return res;
        }
    };

    template<class Type>
    void printVector(std::vector<Type> &v){
        for(const auto &a : v){
            std::cout << a << ", ";
        }
        std::cout << std::endl;
    }


    template<class Type>
    void printSize(std::vector<Type> &v){
        for(const auto &a : v){
            std::cout << a.size() << ", ";
        }
        std::cout << std::endl;
    }



    class h_sum_path_intersection {


    private:
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<std::vector<uint64_t>> m_intersection;
        std::vector<double> m_r;
        std::vector<double> m_l;
        uint64_t m_max_p;
        uint64_t m_sigma;


    public:
        h_sum_path_intersection(const std::vector<PairPredPos> &preds,
                                    bwt_nose &L_S, bwt_type &wt_pred_s,
                                    uint64_t maxP, uint64_t sigma,
                                    uint64_t n_predicates,
                                    query_type q_type = var_var) {


            //auto t0 = std::chrono::high_resolution_clock::now();
            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            for (uint64_t i = 0; i < preds.size(); ++i) {
                const auto& pair = preds[i];
                auto v_target = distinct_values(L_S.get_C(pair.pred_e), L_S.get_C(pair.pred_e + 1) - 1, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.pred_b, m_max_p);
                auto v_source = distinct_values(L_S.get_C(rev_id), L_S.get_C(rev_id + 1) - 1, wt_pred_s);
                m_s.push_back(v_source);

                if(i < preds.size()-1 && pair.pos_e == preds[i+1].pos_b-1){
                    std::vector<std::array<uint64_t, 2ul>> ranges;
                    auto Is_p1 = std::pair<uint64_t, uint64_t>(L_S.get_C(pair.pred_e), L_S.get_C(pair.pred_e + 1) - 1);
                    auto Is_p2 = std::pair<uint64_t, uint64_t>(L_S.get_C(reverse(preds[i+1].pred_b, m_max_p)),
                                                               L_S.get_C(reverse(preds[i+1].pred_b, m_max_p) + 1) - 1);
                    ranges.push_back({Is_p1.first, Is_p1.second});
                    ranges.push_back({Is_p2.first, Is_p2.second});
                    m_intersection.emplace_back(L_S.intersect_nofreq(ranges));
                }else{
                    m_intersection.emplace_back(std::vector<uint64_t>());
                }
            }


            /*auto t1 = std::chrono::high_resolution_clock::now();
            auto intersections = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
            std::cout << "Time intersections: " << intersections << std::endl;*/
            // m_s.push_back(-1ULL);
            auto s = m_s.size();
            std::vector<double> m_d(s), m_i(s);

            for (uint64_t i = 0; i < s; ++i) {
                m_d[i] = m_t[i] / (double) m_s[i];
                m_i[i] = m_s[i] / (double) m_t[i];
            }
            //m_sol_l.resize(s);
            //m_sol_r.resize(s);
            //m_sol_r[s - 1] = m_d[s - 1];
            //m_sol_l[0] = m_i[0];
            m_l.resize(s);
            m_r.resize(s);
            m_r[s - 1] = m_d[s - 1];
            m_l[0] = m_i[0];
            for (uint64_t i = 1; i < s; ++i) {
                //Solutions:
                // m_sol_r[i]: multiplicity ratio of solutions from i to s-1
                // m_sol_l[i]: multiplicity ratio of solutions from 0 to i
                //m_sol_r[s - 1 - i] = m_sol_r[s - i] * m_d[s - 1 - i];
                //m_sol_l[i] = m_sol_l[i - 1] * m_i[i];

                //Paths:
                //m_r: multiplicity ratio of paths from i to s-1
                //m_l: multiplicity ratio of paths from 0 to i
                //m_r[s - 1 - i] = (1 + m_r[s - i]) * m_d[s - 1 - i];
                //m_l[i] = (1 + m_l[i - 1]) * m_i[i];

                m_r[s - 1 - i] = (1 +m_r[s - i]) * m_d[s - 1 - i];
                m_l[i] = (1 + m_l[i - 1]) * m_i[i];
            }

            if(q_type == const_var && preds[0].pos_b == 1){
                m_s[0] = 1;
            }

            if(q_type == var_const && preds[s-1].pos_e == n_predicates){
                m_t[s-1] = 1;
            }

            /*std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);
            std::cout << "----Int----" << std::endl;
            printSize(m_intersection);*/
            /*auto t1 = std::chrono::high_resolution_clock::now();
            auto intersections = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
            std::cout << "Decision: " << intersections << std::endl;*/
            std::cout << "Mandatory" << std::endl;

        }


        info simple(const uint64_t ith) {
            info res;
            double w_left, w_right;
            //std::cout << "Target with ith=" << ith << " size=" << m_t[ith] << std::endl;
            //std::cout << "Source with ith=" << ith << " size=" << m_s[ith] << std::endl;
            if (m_s[ith] < m_t[ith]) {
                res.split = source;
                if (ith == 0) {
                    //Seed * (1+PathsFactorRight)
                    //res.weight = m_s[ith] * (1 + m_r[ith]);
                    res.weight = m_s[ith] * (1 + m_r[ith]);
                    res.first_left = false;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                //first_right = m_s[ith] * ((1 + m_r[ith]) + m_sol_r[ith] * (1 + m_l[ith - 1]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                //first_left = m_s[ith] * ((1 + m_l[ith - 1]) + m_sol_l[ith - 1] * (1 + m_r[ith]));
                //w_right = m_s[ith] * (1 + m_r[ith]);
                //w_left  = m_s[ith] * (1 + m_l[ith - 1]);
                w_right = m_s[ith] * (1 +m_r[ith]);
                w_left  = m_s[ith] * (1 +m_l[ith - 1]);
            } else {
                res.split = target;
                if (ith == m_t.size() - 1) {
                    //Seed * (1+PathsFactorLeft)
                    //res.weight = m_t[ith] * (1 + m_l[ith]);
                    res.weight = m_t[ith] * ( 1 + m_l[ith]);
                    res.first_left = true;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                //first_right = m_t[ith] * ((1 + m_r[ith + 1]) + m_sol_r[ith + 1] * (1 + m_l[ith]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                //first_left = m_t[ith] * ((1 + m_l[ith]) + m_sol_l[ith] * (1 + m_r[ith + 1]));
                //w_right = m_t[ith] * (1 + m_r[ith+1]);
                //w_left  = m_t[ith] * (1 + m_l[ith]);
                w_right = m_t[ith] * (1 + m_r[ith+1]);
                w_left  = m_t[ith] * (1 + m_l[ith]);

            }
            if (w_left <= w_right) {
                res.weight = w_left + w_right;
                res.first_left = true;
            } else {
                res.weight = w_left + w_right;
                res.first_left = false;
            }
            return res;
        }

        info intersection(const uint64_t ith) {
            info res;
            res.split = intersect;
            double seed, w_right, w_left;
            /*if (m_s[ith + 1] < m_t[ith]) {
                //double p = m_t[ith] / (double) (m_sigma);
                seed = m_s[ith + 1]; //Seed
            } else {
                //double p = m_s[ith+1] / (double) (m_sigma);
                seed = m_t[ith]; //Seed
            }*/
            seed = m_intersection[ith].size();
            //std::cout << "Intersection size: " << seed << std::endl;
            //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
            //first_right = seed * ((1 + m_r[ith + 1]) + m_sol_r[ith + 1] * (1 + m_l[ith]));
            //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
            //first_left = seed * ((1 + m_l[ith]) + m_sol_l[ith] * (1 + m_r[ith + 1]));
            //w_right = seed * (1 + m_r[ith+1]);
            w_right = seed * (1 +m_r[ith+1]);
            //w_left  = seed * (1 + m_l[ith]);
            w_left  = seed * (1 +m_l[ith]);
            if (w_left <= w_right) {
                res.weight = w_left + w_right;
                res.first_left = true;
            } else {
                res.weight = w_left + w_right;
                res.first_left = false;
            }
            return res;
        }

        inline std::vector<uint64_t> get_elements_intersection(uint64_t i){
            return m_intersection[i];
        }
    };




}

#endif //RING_RPQ_SELECTIVITY_HPP
