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

    enum split_type {target, source, intersect};

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

                auto e_d = L_S.get_C(pair.id_pred + 1) - 1;
                auto b_d = L_S.get_C(pair.id_pred);
                auto v_target = distinct_values(b_d, e_d, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.id_pred, m_max_p);
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

    class h_distinct_path {


    private:
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<uint64_t> m_r;
        std::vector<uint64_t> m_l;
        uint64_t m_max_p;
        uint64_t m_sigma;


    public:
        h_distinct_path(const std::vector<PairPredPos> &preds,
                        bwt_nose &L_S, bwt_type &wt_pred_s,
                        uint64_t maxP, uint64_t sigma){

            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            for(const auto &pair : preds){

                auto e_d = L_S.get_C(pair.id_pred + 1)-1;
                auto b_d = L_S.get_C(pair.id_pred);
                auto v_target = distinct_values(b_d, e_d, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.id_pred, m_max_p);
                auto e_r = L_S.get_C(rev_id + 1)-1;
                auto b_r = L_S.get_C(rev_id);
                auto v_source = distinct_values(b_r, e_r, wt_pred_s);
                m_s.push_back(v_source);
            }
           // m_s.push_back(-1ULL);
            auto s = m_s.size()+1;
            m_r.resize(s);
            m_r[s-1]=1;
            m_l.resize(s);
            m_l[0]=1;
            for(uint64_t i = 1; i < s; ++i){
                m_r[s-i-1] = m_r[s-i] * m_t[m_t.size()-i];
                m_l[i] = m_l[i-1] * m_s[i-1];
            }
            /*std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);*/


        }

        info simple(const uint64_t ith){
            info res;
            //double a;
            if(m_s[ith] < m_t[ith]){
                res.split = source;
                if(ith == 0){
                    res.weight = (m_s[ith]/ (double) (m_sigma)) * m_r[ith];
                }else{
                    res.weight = (m_s[ith]/ (double) (m_sigma)) * m_r[ith]
                                 + (m_s[ith]/ (double) (m_sigma)) * m_l[ith];
                }
            }else{
                res.split = target;
                if(ith == m_t.size()-1){
                    res.weight = (m_t[ith]/ (double) (m_sigma)) * m_l[ith];
                }else{
                    res.weight = (m_t[ith]/ (double) (m_sigma)) * m_l[ith]
                                 + (m_t[ith]/ (double) (m_sigma)) * m_r[ith];
                }

            }
            return res;
        }

        info intersection(const uint64_t ith) {
            info res;
            res.split = intersect;
            double base = (double) (m_s[ith+1] * m_t[ith]) / (double) (m_sigma * m_sigma);
            res.weight = m_l[ith] * base + base * m_r[ith];
            return res;
        }

    };

    class h_growup_path {


    private:
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<double> m_r;
        std::vector<double> m_l;
        uint64_t m_max_p;
        uint64_t m_sigma;


    public:
        h_growup_path(const std::vector<PairPredPos> &preds,
                      bwt_nose &L_S, bwt_type &wt_pred_s,
                        uint64_t maxP, uint64_t sigma){

            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            for(const auto &pair : preds){

                auto e_d = L_S.get_C(pair.id_pred + 1)-1;
                auto b_d = L_S.get_C(pair.id_pred);
                auto v_target = distinct_values(b_d, e_d, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.id_pred, m_max_p);
                auto e_r = L_S.get_C(rev_id + 1)-1;
                auto b_r = L_S.get_C(rev_id);
                auto v_source = distinct_values(b_r, e_r, wt_pred_s);
                m_s.push_back(v_source);
            }
            // m_s.push_back(-1ULL);
            auto s = m_s.size();
            m_r.resize(s);
            //m_r[s-1]=1;
            m_l.resize(s);
           // m_l[0]=1;
            for(uint64_t i = 0; i < s; ++i){
                m_r[i] = m_t[i] / (double) m_s[i];
                m_l[i] = m_s[i] / (double) m_t[i];
            }
            std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);


        }

        info simple(const uint64_t ith){
            info res;
            //double a;
            double b_l, b_r, w_l, w_r;
            if(m_s[ith] < m_t[ith]){
                res.split = source;
                b_l = b_r = m_s[ith];
                //Right part
                w_r = b_r;
                for(uint64_t i = ith; i < m_r.size(); ++i){
                    b_r = b_r * m_r[i];
                    w_r += b_r;
                }
                //Left part
                w_l = b_l; //Jump from target to source
                for(int64_t i = ith-1; i >= 0; --i){
                    b_l = b_l * m_l[i];
                    w_l += b_l;
                }
            }else{
                res.split = target;
                b_l = b_r = m_t[ith];

                //Right part
                w_r = b_r;
                for(uint64_t i = ith+1; i < m_r.size(); ++i){
                    b_r = b_r * m_r[i];
                    w_r += b_r;
                }
                //Left part
                w_l = b_l; //Jump from target to source
                for(int64_t i = ith; i >= 0; --i){
                    b_l = b_l * m_l[i];
                    w_l += b_l;
                }
            }
            if(w_l > w_r) res.first_left = false;
            res.weight = w_l + w_r;
            return res;
        }

        info intersection(const uint64_t ith) {
            info res;
            res.split = intersect;
            double b_l, b_r, w_l, w_r;
            if(m_s[ith+1] < m_t[ith]){
                double p = m_t[ith] / (double) (m_sigma);
                b_l = b_r = p * m_s[ith+1];
            }else{
                double p = m_s[ith+1] / (double) (m_sigma);
                b_l = b_r = p * m_t[ith];
            }
            //Right part
            w_r = b_r; //Jump from source to target
            for(uint64_t i = ith+1; i < m_r.size(); ++i){
                b_r = b_r * m_r[i];
                w_r += b_r;
            }
            //Left part
            w_l = b_l; //Jump from target to source
            for(int64_t i = ith; i >= 0; --i){
                b_l = b_l * m_l[i];
                w_l += b_l;
            }
            if(w_l > w_r) res.first_left = false;
            res.weight = w_l + w_r;
            return res;
        }

    };


    class h_growup_path2 {


    private:
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<double> m_r;
        std::vector<double> m_l;
        uint64_t m_max_p;
        uint64_t m_sigma;


    public:
        h_growup_path2(const std::vector<PairPredPos> &preds,
                       bwt_nose &L_S, bwt_type &wt_pred_s,
                      uint64_t maxP, uint64_t sigma){

            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            for(const auto &pair : preds){

                auto e_d = L_S.get_C(pair.id_pred + 1)-1;
                auto b_d = L_S.get_C(pair.id_pred);
                auto v_target = distinct_values(b_d, e_d, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.id_pred, m_max_p);
                auto e_r = L_S.get_C(rev_id + 1)-1;
                auto b_r = L_S.get_C(rev_id);
                auto v_source = distinct_values(b_r, e_r, wt_pred_s);
                m_s.push_back(v_source);
            }
            // m_s.push_back(-1ULL);
            auto s = m_s.size();
            m_r.resize(s);
            //m_r[s-1]=1;
            m_l.resize(s);
            // m_l[0]=1;
            for(uint64_t i = 0; i < s; ++i){
                m_r[i] = m_t[i] / (double) m_s[i];
                m_l[i] = m_s[i] / (double) m_t[i];
            }
            /*std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);*/


        }

        info simple(const uint64_t ith){
            info res;
            //double a;
            double b1_l, b1_r, b2_l, b2_r, w1_l, w2_l, w1_r, w2_r;
            if(m_s[ith] < m_t[ith]){
                res.split = source;
                b1_l = b1_r = m_s[ith];
                //Right part as first option
                w1_r = b1_r;
                for(uint64_t i = ith; i < m_r.size(); ++i){
                    b1_r = b1_r * m_r[i];
                    w1_r += b1_r;
                }
                //Left part as first option
                w1_l = b1_l; //Jump from target to source
                for(int64_t i = ith-1; i >= 0; --i){
                    b1_l = b1_l * m_l[i];
                    w1_l += b1_l;
                }
                //Left part as second option
                b2_l = b1_r;
                w2_l = b2_l; //Jump from target to source
                for(int64_t i = ith-1; i >= 0; --i){
                    b2_l = b2_l * m_l[i];
                    w2_l += b2_l;
                }
                //Right part as second option
                b2_r = b1_l;
                w2_r = b2_r;
                for(uint64_t i = ith; i < m_r.size(); ++i){
                    b2_r = b2_r * m_r[i];
                    w2_r += b2_r;
                }

            }else{
                res.split = target;
                b1_l = b1_r = m_t[ith];
                //Right part as first option
                w1_r = b1_r;
                for(uint64_t i = ith+1; i < m_r.size(); ++i){
                    b1_r = b1_r * m_r[i];
                    w1_r += b1_r;
                }
                //Left part as first option
                w1_l = b1_l; //Jump from target to source
                for(int64_t i = ith; i >= 0; --i){
                    b1_l = b1_l * m_l[i];
                    w1_l += b1_l;
                }
                //Left part as second option
                b2_l = b1_r;
                w2_l = b2_l; //Jump from target to source
                for(int64_t i = ith; i >= 0; --i){
                    b2_l = b2_l * m_l[i];
                    w2_l += b2_l;
                }
                //Right part as second option
                b2_r = b1_l;
                w2_r = b2_r;
                for(uint64_t i = ith+1; i < m_r.size(); ++i){
                    b2_r = b2_r * m_r[i];
                    w2_r += b2_r;
                }
            }
            auto first_left = w1_l + w2_r;
            auto first_right = w1_r + w2_l;
            if(first_left <= first_right){
                res.weight = first_left;
                res.first_left = true;
            }else{
                res.weight = first_right;
                res.first_left = false;
            }
            return res;
        }

        info intersection(const uint64_t ith) {
            info res;
            res.split = intersect;
            double b1_l, b1_r, b2_l, b2_r, w1_l, w1_r, w2_l, w2_r;
            if(m_s[ith+1] < m_t[ith]){
                //double p = m_t[ith] / (double) (m_sigma);
                b1_l = b1_r = m_s[ith+1];
            }else{
                //double p = m_s[ith+1] / (double) (m_sigma);
                b1_l = b1_r = m_t[ith];
            }
            //Right part as first option
            w1_r = b1_r; //Jump from source to target
            for(uint64_t i = ith+1; i < m_r.size(); ++i){
                b1_r = b1_r * m_r[i];
                w1_r += b1_r;
            }
            //Left part as first option
            w1_l = b1_l; //Jump from target to source
            for(int64_t i = ith; i >= 0; --i){
                b1_l = b1_l * m_l[i];
                w1_l += b1_l;
            }
            //Left part as second option
            b2_l = b1_r;
            w2_l = b2_l; //Jump from target to source
            for(int64_t i = ith; i >= 0; --i){
                b2_l = b2_l * m_l[i];
                w2_l += b2_l;
            }
            //Right part as second option
            b2_r = b1_l;
            w2_r = b2_r; //Jump from source to target
            for(uint64_t i = ith+1; i < m_r.size(); ++i){
                b2_r = b2_r * m_r[i];
                w2_r += b2_r;
            }
            auto first_left = w1_l + w2_r;
            auto first_right = w1_r + w2_l;
            if(first_left <= first_right){
                res.weight = first_left;
                res.first_left = true;
            }else{
                res.weight = first_right;
                res.first_left = false;
            }
            return res;
        }

    };

    class h_growup_path2_opt {


    private:
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<double> m_r;
        std::vector<double> m_sol_r;
        std::vector<double> m_l;
        std::vector<double> m_sol_l;
        uint64_t m_max_p;
        uint64_t m_sigma;


    public:
        h_growup_path2_opt(const std::vector<PairPredPos> &preds,
                        bwt_nose &L_S, bwt_type &wt_pred_s,
                       uint64_t maxP, uint64_t sigma){

            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            for(const auto &pair : preds){

                auto e_d = L_S.get_C(pair.id_pred + 1)-1;
                auto b_d = L_S.get_C(pair.id_pred);
                auto v_target = distinct_values(b_d, e_d, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.id_pred, m_max_p);
                auto e_r = L_S.get_C(rev_id + 1)-1;
                auto b_r = L_S.get_C(rev_id);
                auto v_source = distinct_values(b_r, e_r, wt_pred_s);
                m_s.push_back(v_source);
            }
            // m_s.push_back(-1ULL);
            auto s = m_s.size();
            std::vector<double> m_d(s), m_i(s);

            for(uint64_t i = 0; i < s; ++i){
                m_d[i] = m_t[i] / (double) m_s[i];
                m_i[i] = m_s[i] / (double) m_t[i];
            }
            m_sol_l.resize(s);
            m_sol_r.resize(s);
            m_sol_r[s-1] = m_d[s-1];
            m_sol_l[0] = m_i[0];
            m_l.resize(s);
            m_r.resize(s);
            m_r[s-1] = m_d[s-1];
            m_l[0] = m_i[0];
            for(uint64_t i = 1; i < s; ++i){
                //Solutions:
                // m_sol_r[i]: multiplicity ratio of solutions from i to s-1
                // m_sol_l[i]: multiplicity ratio of solutions from 0 to i
                m_sol_r[s-1-i] = m_sol_r[s-i] * m_d[s-1-i];
                m_sol_l[i] = m_sol_l[i-1] * m_i[i];

                //Paths:
                //m_r: multiplicity ratio of paths from i to s-1
                //m_l: multiplicity ratio of paths from 0 to i
                m_r[s-1-i] = (1+m_r[s-i])*m_d[s-1-i];
                m_l[i] = (1+m_l[i-1])*m_i[i];
            }

           /* std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);*/

        }


        info simple(const uint64_t ith){
            info res;
            double first_left, first_right;
            if(m_s[ith] < m_t[ith]){
                res.split = source;
                if(ith == 0){
                    //Seed * (1+PathsFactorRight)
                    res.weight = m_s[ith] * (1+m_r[ith]);
                    res.first_left = false;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                first_right = m_s[ith] * ((1+m_r[ith]) + m_sol_r[ith] * (1+m_l[ith-1]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                first_left = m_s[ith] * ((1+m_l[ith-1]) + m_sol_l[ith-1] * (1+m_r[ith]));
            }else{
                res.split = target;
                if(ith == m_t.size()-1){
                    //Seed * (1+PathsFactorLeft)
                    res.weight = m_t[ith] * (1+m_l[ith]);
                    res.first_left = true;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                first_right = m_t[ith] * ((1+m_r[ith+1]) + m_sol_r[ith+1] * (1+m_l[ith]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                first_left = m_t[ith] * ((1+m_l[ith]) + m_sol_l[ith] * (1+m_r[ith+1]));
            }
            if(first_left <= first_right){
                res.weight = first_left;
                res.first_left = true;
            }else{
                res.weight = first_right;
                res.first_left = false;
            }
            return res;
        }

        info intersection(const uint64_t ith) {
            info res;
            res.split = intersect;
            double seed, first_right, first_left;
            if(m_s[ith+1] < m_t[ith]){
                //double p = m_t[ith] / (double) (m_sigma);
                seed = m_s[ith+1]; //Seed
            }else{
                //double p = m_s[ith+1] / (double) (m_sigma);
                seed = m_t[ith]; //Seed
            }
            //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
            first_right = seed * ((1+m_r[ith+1]) + m_sol_r[ith+1] * (1+m_l[ith]));
            //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
            first_left = seed * ((1+m_l[ith]) + m_sol_l[ith] * (1+m_r[ith+1]));
            if(first_left <= first_right){
                res.weight = first_left;
                res.first_left = true;
            }else{
                res.weight = first_right;
                res.first_left = false;
            }
            return res;
        }

    };

    class h_growup_path3 {


    private:
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<double> m_r;
        std::vector<double> m_sol_r;
        std::vector<double> m_l;
        std::vector<double> m_sol_l;
        uint64_t m_max_p;
        uint64_t m_sigma;


    public:
        h_growup_path3(const std::vector<PairPredPos> &preds,
                       bwt_nose &L_S, bwt_type &wt_pred_s,
                           uint64_t maxP, uint64_t sigma){

            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            std::vector<double> lengths;
            for(const auto &pair : preds){

                auto e_d = L_S.get_C(pair.id_pred + 1)-1;
                auto b_d = L_S.get_C(pair.id_pred);
                auto v_target = distinct_values(b_d, e_d, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.id_pred, m_max_p);
                auto e_r = L_S.get_C(rev_id + 1)-1;
                auto b_r = L_S.get_C(rev_id);
                auto v_source = distinct_values(b_r, e_r, wt_pred_s);
                m_s.push_back(v_source);

                lengths.push_back(e_r-b_r+1);
            }
            // m_s.push_back(-1ULL);
            auto s = m_s.size();
            std::vector<double> m_d(s), m_i(s);
            for(uint64_t i = 0; i < s; ++i){
                m_d[i] = m_t[i] * (1 - std::pow(1 - (1/(double) m_t[i]), lengths[i]/ (double) m_s[i]));
                m_i[i] = m_s[i] * (1 - std::pow(1 - (1/(double) m_s[i]), lengths[i]/ (double) m_t[i]));
            }
            m_sol_l.resize(s);
            m_sol_r.resize(s);
            m_sol_r[s-1] = m_d[s-1];
            m_sol_l[0] = m_i[0];
            m_l.resize(s);
            m_r.resize(s);
            m_r[s-1] = m_d[s-1];
            m_l[0] = m_i[0];
            for(uint64_t i = 1; i < s; ++i){
                //Solutions:
                // m_sol_r[i]: multiplicity ratio of solutions from i to s-1
                // m_sol_l[i]: multiplicity ratio of solutions from 0 to i
                m_sol_r[s-1-i] = m_sol_r[s-i] * m_d[s-1-i];
                m_sol_l[i] = m_sol_l[i-1] * m_i[i];

                //Paths:
                //m_r: multiplicity ratio of paths from i to s-1
                //m_l: multiplicity ratio of paths from 0 to i
                m_r[s-1-i] = (1+m_r[s-i])*m_d[s-1-i];
                m_l[i] = (1+m_l[i-1])*m_i[i];
            }

            /*std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);*/

        }

        info simple(const uint64_t ith){
            info res;
            double first_left, first_right;
            if(m_s[ith] < m_t[ith]){
                res.split = source;
                if(ith == 0){
                    //Seed * (1+PathsFactorRight)
                    res.weight = m_s[ith] * (1+m_r[ith]);
                    res.first_left = false;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                first_right = m_s[ith] * ((1+m_r[ith]) + m_sol_r[ith] * (1+m_l[ith-1]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                first_left = m_s[ith] * ((1+m_l[ith-1]) + m_sol_l[ith-1] * (1+m_r[ith]));
            }else{
                res.split = target;
                if(ith == m_t.size()-1){
                    //Seed * (1+PathsFactorLeft)
                    res.weight = m_t[ith] * (1+m_l[ith]);
                    res.first_left = true;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                first_right = m_t[ith] * ((1+m_r[ith+1]) + m_sol_r[ith+1] * (1+m_l[ith]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                first_left = m_t[ith] * ((1+m_l[ith]) + m_sol_l[ith] * (1+m_r[ith+1]));
            }
            if(first_left <= first_right){
                res.weight = first_left;
                res.first_left = true;
            }else{
                res.weight = first_right;
                res.first_left = false;
            }
            return res;
        }

        info intersection(const uint64_t ith) {
            info res;
            res.split = intersect;
            double seed, first_right, first_left;
            if(m_s[ith+1] < m_t[ith]){
                //double p = m_t[ith] / (double) (m_sigma);
                seed = m_s[ith+1]; //Seed
            }else{
                //double p = m_s[ith+1] / (double) (m_sigma);
                seed = m_t[ith]; //Seed
            }
            //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
            first_right = seed * ((1+m_r[ith+1]) + m_sol_r[ith+1] * (1+m_l[ith]));
            //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
            first_left = seed * ((1+m_l[ith]) + m_sol_l[ith] * (1+m_r[ith+1]));
            if(first_left <= first_right){
                res.weight = first_left;
                res.first_left = true;
            }else{
                res.weight = first_right;
                res.first_left = false;
            }
            return res;
        }

    };


    class h_sum_path2_intersection {


    private:
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<std::vector<uint64_t>> m_intersection;
        std::vector<double> m_r;
        std::vector<double> m_sol_r;
        std::vector<double> m_l;
        std::vector<double> m_sol_l;
        uint64_t m_max_p;
        uint64_t m_sigma;


    public:
        h_sum_path2_intersection(const std::vector<PairPredPos> &preds,
                                    bwt_nose &L_S, bwt_type &wt_pred_s,
                           uint64_t maxP, uint64_t sigma) {


            auto t0 = std::chrono::high_resolution_clock::now();
            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            for (uint64_t i = 0; i < preds.size(); ++i) {
                const auto& pair = preds[i];
                auto v_target = distinct_values(L_S.get_C(pair.id_pred), L_S.get_C(pair.id_pred + 1) - 1, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.id_pred, m_max_p);
                auto v_source = distinct_values(L_S.get_C(rev_id), L_S.get_C(rev_id + 1) - 1, wt_pred_s);
                m_s.push_back(v_source);

                if(i < preds.size()-1 && pair.pos == preds[i+1].pos-1){
                    std::vector<std::array<uint64_t, 2ul>> ranges;
                    auto Is_p1 = std::pair<uint64_t, uint64_t>(L_S.get_C(pair.id_pred), L_S.get_C(pair.id_pred + 1) - 1);
                    auto Is_p2 = std::pair<uint64_t, uint64_t>(L_S.get_C(reverse(preds[i+1].id_pred, m_max_p)),
                                                               L_S.get_C(reverse(preds[i+1].id_pred, m_max_p) + 1) - 1);
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
            m_sol_l.resize(s);
            m_sol_r.resize(s);
            m_sol_r[s - 1] = m_d[s - 1];
            m_sol_l[0] = m_i[0];
            m_l.resize(s);
            m_r.resize(s);
            m_r[s - 1] = m_d[s - 1];
            m_l[0] = m_i[0];
            for (uint64_t i = 1; i < s; ++i) {
                //Solutions:
                // m_sol_r[i]: multiplicity ratio of solutions from i to s-1
                // m_sol_l[i]: multiplicity ratio of solutions from 0 to i
                m_sol_r[s - 1 - i] = m_sol_r[s - i] * m_d[s - 1 - i];
                m_sol_l[i] = m_sol_l[i - 1] * m_i[i];

                //Paths:
                //m_r: multiplicity ratio of paths from i to s-1
                //m_l: multiplicity ratio of paths from 0 to i
                m_r[s - 1 - i] = (1 + m_r[s - i]) * m_d[s - 1 - i];
                m_l[i] = (1 + m_l[i - 1]) * m_i[i];
            }
            /*
            std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);*/

            auto t1 = std::chrono::high_resolution_clock::now();
            auto intersections = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
            std::cout << "Decision: " << intersections << std::endl;

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
                    res.weight = m_s[ith] * (1 + m_r[ith]);
                    res.first_left = false;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                //first_right = m_s[ith] * ((1 + m_r[ith]) + m_sol_r[ith] * (1 + m_l[ith - 1]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                //first_left = m_s[ith] * ((1 + m_l[ith - 1]) + m_sol_l[ith - 1] * (1 + m_r[ith]));
                w_right = m_s[ith] * (1 + m_r[ith]);
                w_left  = m_s[ith] * (1 + m_l[ith - 1]);
            } else {
                res.split = target;
                if (ith == m_t.size() - 1) {
                    //Seed * (1+PathsFactorLeft)
                    res.weight = m_t[ith] * (1 + m_l[ith]);
                    res.first_left = true;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                //first_right = m_t[ith] * ((1 + m_r[ith + 1]) + m_sol_r[ith + 1] * (1 + m_l[ith]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                //first_left = m_t[ith] * ((1 + m_l[ith]) + m_sol_l[ith] * (1 + m_r[ith + 1]));
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
            w_right = seed * (1 + m_r[ith+1]);
            w_left  = seed * (1 + m_l[ith]);
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


    class h_sum_path3_intersection {


    private:
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<std::vector<uint64_t>> m_intersection;
        std::vector<double> m_r;
        std::vector<double> m_sol_r;
        std::vector<double> m_l;
        std::vector<double> m_sol_l;
        uint64_t m_max_p;
        uint64_t m_sigma;


    public:
        h_sum_path3_intersection(const std::vector<PairPredPos> &preds,
                                 bwt_nose &L_S, bwt_type &wt_pred_s,
                                 uint64_t maxP, uint64_t sigma) {


            auto t0 = std::chrono::high_resolution_clock::now();
            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            for (uint64_t i = 0; i < preds.size(); ++i) {
                const auto& pair = preds[i];
                auto v_target = distinct_values(L_S.get_C(pair.id_pred), L_S.get_C(pair.id_pred + 1) - 1, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.id_pred, m_max_p);
                auto v_source = distinct_values(L_S.get_C(rev_id), L_S.get_C(rev_id + 1) - 1, wt_pred_s);
                m_s.push_back(v_source);

                if(i < preds.size()-1 && pair.pos == preds[i+1].pos-1){
                    std::vector<std::array<uint64_t, 2ul>> ranges;
                    auto Is_p1 = std::pair<uint64_t, uint64_t>(L_S.get_C(pair.id_pred), L_S.get_C(pair.id_pred + 1) - 1);
                    auto Is_p2 = std::pair<uint64_t, uint64_t>(L_S.get_C(reverse(preds[i+1].id_pred, m_max_p)),
                                                               L_S.get_C(reverse(preds[i+1].id_pred, m_max_p) + 1) - 1);
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
            m_sol_l.resize(s);
            m_sol_r.resize(s);
            m_sol_r[s - 1] = m_d[s - 1];
            m_sol_l[0] = m_i[0];
            m_l.resize(s);
            m_r.resize(s);
            m_r[s - 1] = m_d[s - 1];
            m_l[0] = m_i[0];
            for (uint64_t i = 1; i < s; ++i) {
                //Solutions:
                // m_sol_r[i]: multiplicity ratio of solutions from i to s-1
                // m_sol_l[i]: multiplicity ratio of solutions from 0 to i
                m_sol_r[s - 1 - i] = m_sol_r[s - i] * m_d[s - 1 - i];
                m_sol_l[i] = m_sol_l[i - 1] * m_i[i];

                //Paths:
                //m_r: multiplicity ratio of paths from i to s-1
                //m_l: multiplicity ratio of paths from 0 to i
                m_r[s - 1 - i] = (1 + m_r[s - i]) * m_d[s - 1 - i];
                m_l[i] = (1 + m_l[i - 1]) * m_i[i];
            }
            /*
            std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);*/

            auto t1 = std::chrono::high_resolution_clock::now();
            auto intersections = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
            std::cout << "Decision: " << intersections << std::endl;

        }


        info simple(const uint64_t ith) {
            info res;
            double first_left, first_right;
            //std::cout << "Target with ith=" << ith << " size=" << m_t[ith] << std::endl;
            //std::cout << "Source with ith=" << ith << " size=" << m_s[ith] << std::endl;
            if (m_s[ith] < m_t[ith]) {
                res.split = source;
                if (ith == 0) {
                    //Seed * (1+PathsFactorRight)
                    res.weight = m_s[ith] * (1 + m_r[ith]);
                    res.first_left = false;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                double seed = m_s[ith];
                first_right = seed * ((1 + m_r[ith]) + std::min(seed, m_sol_r[ith]) * (1 + m_l[ith - 1]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                first_left = seed * ((1 + m_l[ith - 1]) + std::min(seed, m_sol_l[ith - 1]) * (1 + m_r[ith]));

            } else {
                res.split = target;
                if (ith == m_t.size() - 1) {
                    //Seed * (1+PathsFactorLeft)
                    res.weight = m_t[ith] * (1 + m_l[ith]);
                    res.first_left = true;
                    return res;
                }
                double seed = m_t[ith];
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                first_right = seed * ((1 + m_r[ith + 1]) + std::min(seed, m_sol_r[ith + 1]) * (1 + m_l[ith]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                first_left = seed * ((1 + m_l[ith]) + std::min(seed, m_sol_l[ith]) * (1 + m_r[ith + 1]));

            }
            if (first_left <= first_right) {
                res.weight = first_left;
                res.first_left = true;
            } else {
                res.weight = first_right;
                res.first_left = false;
            }
            return res;
        }

        info intersection(const uint64_t ith) {
            info res;
            res.split = intersect;
            double seed, first_right, first_left;
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
            first_right = seed * ((1 + m_r[ith + 1]) + std::min(seed, m_sol_r[ith + 1]) * (1 + m_l[ith]));
            //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
            first_left = seed * ((1 + m_l[ith - 1]) + std::min(seed, m_sol_l[ith - 1]) * (1 + m_r[ith]));
            if (first_left <= first_right) {
                res.weight = first_left;
                res.first_left = true;
            } else {
                res.weight = first_right;
                res.first_left = false;
            }
            return res;
        }

        inline std::vector<uint64_t> get_elements_intersection(uint64_t i){
            return m_intersection[i];
        }
    };

    class h_sum_path_intersection {


    private:
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<std::vector<uint64_t>> m_intersection;
        std::vector<double> m_r;
        std::vector<double> m_sol_r;
        std::vector<double> m_l;
        std::vector<double> m_sol_l;
        uint64_t m_max_p;
        uint64_t m_sigma;
        const std::vector<PairPredPos>* m_preds;


    public:
        h_sum_path_intersection(const std::vector<PairPredPos> &preds,
                                 bwt_nose &L_S, bwt_type &wt_pred_s,
                                 uint64_t maxP, uint64_t sigma) {

            auto t0 = std::chrono::high_resolution_clock::now();
            m_preds = &preds;
            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            for (uint64_t i = 0; i < preds.size(); ++i) {
                const auto& pair = preds[i];
                auto v_target = distinct_values(L_S.get_C(pair.id_pred), L_S.get_C(pair.id_pred + 1) - 1, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.id_pred, m_max_p);
                auto v_source = distinct_values(L_S.get_C(rev_id), L_S.get_C(rev_id + 1) - 1, wt_pred_s);
                m_s.push_back(v_source);

                if(i < preds.size()-1 && pair.pos == preds[i+1].pos-1){
                    std::vector<std::array<uint64_t, 2ul>> ranges;
                    auto Is_p1 = std::pair<uint64_t, uint64_t>(L_S.get_C(pair.id_pred), L_S.get_C(pair.id_pred + 1) - 1);
                    auto Is_p2 = std::pair<uint64_t, uint64_t>(L_S.get_C(reverse(preds[i+1].id_pred, m_max_p)),
                                                               L_S.get_C(reverse(preds[i+1].id_pred, m_max_p) + 1) - 1);
                    ranges.push_back({Is_p1.first, Is_p1.second});
                    ranges.push_back({Is_p2.first, Is_p2.second});
                    m_intersection.emplace_back(L_S.intersect_nofreq(ranges));
                }else{
                    m_intersection.emplace_back(std::vector<uint64_t>());
                }
            }
            // m_s.push_back(-1ULL);
            auto s = m_s.size();
            std::vector<double> m_d(s), m_i(s);

            for (uint64_t i = 0; i < s; ++i) {
                m_d[i] = m_t[i] / (double) m_s[i];
                m_i[i] = m_s[i] / (double) m_t[i];
            }
            m_sol_l.resize(s);
            m_sol_r.resize(s);
            m_sol_r[s - 1] = m_d[s - 1];
            m_sol_l[0] = m_i[0];
            m_l.resize(s);
            m_r.resize(s);
            m_r[s - 1] = m_d[s - 1];
            m_l[0] = m_i[0];
            for (uint64_t i = 1; i < s; ++i) {
                //Solutions:
                // m_sol_r[i]: multiplicity ratio of solutions from i to s-1
                // m_sol_l[i]: multiplicity ratio of solutions from 0 to i
                m_sol_r[s - 1 - i] = m_sol_r[s - i] * m_d[s - 1 - i];
                m_sol_l[i] = m_sol_l[i - 1] * m_i[i];

                //Paths:
                //m_r: multiplicity ratio of paths from i to s-1
                //m_l: multiplicity ratio of paths from 0 to i
                m_r[s - 1 - i] = (1 + m_r[s - i]) * m_d[s - 1 - i];
                m_l[i] = (1 + m_l[i - 1]) * m_i[i];
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            auto intersections = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
            std::cout << "Decision: " << intersections << std::endl;
            /*
            std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);*/

        }


        info_preds simple(const uint64_t ith) {
            info_preds res;
            double first_left, first_right;
            if (m_s[ith] < m_t[ith]) {
                res.split = source;
                if (ith == 0) {
                    //Seed * (1+PathsFactorRight)
                    res.weight = m_s[ith] * (1 + m_r[ith]);
                    res.mand_pred_left = m_max_p+1;
                    res.mand_pred_right = m_preds->at(ith).id_pred;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                first_right = m_s[ith] * ((1 + m_r[ith]) + m_sol_r[ith] * (1 + m_l[ith - 1]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                first_left = m_s[ith] * ((1 + m_l[ith - 1]) + m_sol_l[ith - 1] * (1 + m_r[ith]));
                res.mand_pred_left =  reverse(m_preds->at(ith-1).id_pred, m_max_p);
                res.mand_pred_right = m_preds->at(ith).id_pred;
            } else {
                res.split = target;
                if (ith == m_t.size() - 1) {
                    //Seed * (1+PathsFactorLeft)
                    res.weight = m_t[ith] * (1 + m_l[ith]);
                    res.mand_pred_left =  reverse(m_preds->at(ith).id_pred, m_max_p);
                    res.mand_pred_right = m_max_p+1;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                first_right = m_t[ith] * ((1 + m_r[ith + 1]) + m_sol_r[ith + 1] * (1 + m_l[ith]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                first_left = m_t[ith] * ((1 + m_l[ith]) + m_sol_l[ith] * (1 + m_r[ith + 1]));
                res.mand_pred_left =  reverse(m_preds->at(ith).id_pred, m_max_p);
                res.mand_pred_right = m_preds->at(ith+1).id_pred;
            }
            if (first_left <= first_right) {
                res.weight = first_left;
            } else {
                res.weight = first_right;
            }
            return res;
        }

        info_preds intersection(const uint64_t ith) {
            info_preds res;
            res.split = intersect;
            double seed, first_right, first_left;
            /*if (m_s[ith + 1] < m_t[ith]) {
                //double p = m_t[ith] / (double) (m_sigma);
                seed = m_s[ith + 1]; //Seed
            } else {
                //double p = m_s[ith+1] / (double) (m_sigma);
                seed = m_t[ith]; //Seed
            }*/
            seed = m_intersection[ith].size();
            //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
            first_right = seed * ((1 + m_r[ith + 1]) + m_sol_r[ith + 1] * (1 + m_l[ith]));
            //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
            first_left = seed * ((1 + m_l[ith]) + m_sol_l[ith] * (1 + m_r[ith + 1]));
            if (first_left <= first_right) {
                res.weight = first_left;
            } else {
                res.weight = first_right;
            }
            res.mand_pred_left =  reverse(m_preds->at(ith).id_pred, m_max_p);
            res.mand_pred_right = m_preds->at(ith+1).id_pred;
            return res;
        }

        inline std::vector<uint64_t> get_elements_intersection(uint64_t i){
            return m_intersection[i];
        }
    };

    class h_distinct_intersection {


    private:
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<std::vector<uint64_t>> m_intersection;
        uint64_t m_max_p;
        uint64_t m_sigma;
        const std::vector<PairPredPos>* m_preds;


    public:
        h_distinct_intersection(const std::vector<PairPredPos> &preds,
                                bwt_nose &L_S, bwt_type &wt_pred_s,
                                uint64_t maxP, uint64_t sigma) {

            m_preds = &preds;
            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            auto t0 = std::chrono::high_resolution_clock::now();
            for (const auto& pair : preds) {
                auto v_target = distinct_values(L_S.get_C(pair.id_pred),
                                                L_S.get_C(pair.id_pred + 1) - 1, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.id_pred, m_max_p);
                auto v_source = distinct_values(L_S.get_C(rev_id),
                                                L_S.get_C(rev_id + 1) - 1, wt_pred_s);
                m_s.push_back(v_source);

            }

            auto t1 = std::chrono::high_resolution_clock::now();
            for (uint64_t i = 0; i < preds.size(); ++i) {
                if (i < preds.size() - 1 && preds[i].pos == preds[i + 1].pos - 1) {
                    std::vector<std::array<uint64_t, 2ul>> ranges;
                    auto Is_p1 = std::pair<uint64_t, uint64_t>(L_S.get_C(preds[i].id_pred),
                                                               L_S.get_C(preds[i].id_pred + 1) - 1);
                    auto Is_p2 = std::pair<uint64_t, uint64_t>(L_S.get_C(reverse(preds[i + 1].id_pred, m_max_p)),
                                                               L_S.get_C(reverse(preds[i + 1].id_pred, m_max_p) + 1) -
                                                               1);
                    ranges.push_back({Is_p1.first, Is_p1.second});
                    ranges.push_back({Is_p2.first, Is_p2.second});
                    m_intersection.emplace_back(L_S.intersect_nofreq(ranges));
                } else {
                    m_intersection.emplace_back(std::vector<uint64_t>());
                }
            }

            auto t2 = std::chrono::high_resolution_clock::now();
            auto distincts = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
            auto intersections = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
            auto total = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t0).count();
            std::cout << "Distincts: " << distincts << std::endl;
            std::cout << "Intersections: " << intersections << std::endl;
            std::cout << "Decision: " << total << std::endl;
            /*
            std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);*/

        }


        info_preds simple(const uint64_t ith) {
            info_preds res;
            double first_left, first_right;
            if (m_s[ith] < m_t[ith]) {
                res.split = source;
                res.weight = m_s[ith];
                if (ith == 0) {
                    //Seed * (1+PathsFactorRight)
                    res.mand_pred_left = m_max_p+1;
                    res.mand_pred_right = m_preds->at(ith).id_pred;
                    return res;
                }
                res.mand_pred_left =  reverse(m_preds->at(ith-1).id_pred, m_max_p);
                res.mand_pred_right = m_preds->at(ith).id_pred;
            } else {
                res.split = target;
                res.weight = m_t[ith];
                if (ith == m_t.size() - 1) {
                    //Seed * (1+PathsFactorLeft)
                    res.mand_pred_left =  reverse(m_preds->at(ith).id_pred, m_max_p);
                    res.mand_pred_right = m_max_p+1;
                    return res;
                }
                res.mand_pred_left =  reverse(m_preds->at(ith-1).id_pred, m_max_p);
                res.mand_pred_right = m_preds->at(ith).id_pred;
            }
            return res;
        }

        info_preds intersection(const uint64_t ith) {
            info_preds res;
            res.split = intersect;
            res.weight = m_intersection[ith].size();
            res.mand_pred_left =  reverse(m_preds->at(ith).id_pred, m_max_p);
            res.mand_pred_right = m_preds->at(ith+1).id_pred;
            return res;
        }

        inline std::vector<uint64_t> get_elements_intersection(uint64_t i){
            return m_intersection[i];
        }
    };

/*
    class h_sum_path3_intersection {


    private:
        std::vector<uint64_t> m_s;
        std::vector<uint64_t> m_t;
        std::vector<uint64_t> m_intersection;
        std::vector<double> m_r;
        std::vector<double> m_sol_r;
        std::vector<double> m_l;
        std::vector<double> m_sol_l;
        uint64_t m_max_p;
        uint64_t m_sigma;


    public:
        h_sum_path3_intersection(const std::vector<PairPredPos> &preds,
                                    bwt_nose &L_S, bwt_type &wt_pred_s,
                                    uint64_t maxP, uint64_t sigma){

            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            std::vector<double> lengths;
            for (uint64_t i = 0; i < preds.size(); ++i) {
                const auto& pair = preds[i];
                auto e_d = L_S.get_C(pair.id_pred + 1) - 1;
                auto b_d = L_S.get_C(pair.id_pred);
                auto v_target = distinct_values(b_d, e_d, wt_pred_s);
                m_t.push_back(v_target);

                auto rev_id = reverse(pair.id_pred, m_max_p);
                auto e_r = L_S.get_C(rev_id + 1) - 1;
                auto b_r = L_S.get_C(rev_id);
                auto v_source = distinct_values(b_r, e_r, wt_pred_s);
                m_s.push_back(v_source);

                lengths.push_back(e_r-b_r+1);

                if(i < preds.size()-1 && pair.pos == preds[i+1].pos-1){
                    std::vector<std::array<uint64_t, 2ul>> ranges;
                    auto Is_p1 = std::pair<uint64_t, uint64_t>(b_d, e_d);
                    auto Is_p2 = std::pair<uint64_t, uint64_t>(L_S.get_C(reverse(preds[i+1].id_pred, m_max_p)),
                                                               L_S.get_C(reverse(preds[i+1].id_pred, m_max_p) + 1) - 1);
                    ranges.push_back({Is_p1.first, Is_p1.second});
                    ranges.push_back({Is_p2.first, Is_p2.second});
                    m_intersection.push_back(L_S.intersect(ranges).size());
                }else{
                    m_intersection.push_back(0);
                }
            }

            //auto t1 = std::chrono::high_resolution_clock::now();
            //auto intersections = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
            //std::cout << "Time intersections: " << intersections << std::endl;

            // m_s.push_back(-1ULL);
            auto s = m_s.size();
            std::vector<double> m_d(s), m_i(s);
            for(uint64_t i = 0; i < s; ++i){
                m_d[i] = m_t[i] * (1 - std::pow(1 - (1/(double) m_t[i]), lengths[i]/ (double) m_s[i]));
                m_i[i] = m_s[i] * (1 - std::pow(1 - (1/(double) m_s[i]), lengths[i]/ (double) m_t[i]));
            }
            m_sol_l.resize(s);
            m_sol_r.resize(s);
            m_sol_r[s-1] = m_d[s-1];
            m_sol_l[0] = m_i[0];
            m_l.resize(s);
            m_r.resize(s);
            m_r[s-1] = m_d[s-1];
            m_l[0] = m_i[0];
            for(uint64_t i = 1; i < s; ++i){
                //Solutions:
                // m_sol_r[i]: multiplicity ratio of solutions from i to s-1
                // m_sol_l[i]: multiplicity ratio of solutions from 0 to i
                m_sol_r[s-1-i] = m_sol_r[s-i] * m_d[s-1-i];
                m_sol_l[i] = m_sol_l[i-1] * m_i[i];

                //Paths:
                //m_r: multiplicity ratio of paths from i to s-1
                //m_l: multiplicity ratio of paths from 0 to i
                m_r[s-1-i] = (1+m_r[s-i])*m_d[s-1-i];
                m_l[i] = (1+m_l[i-1])*m_i[i];
            }
            
            std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);

        }

        info simple(const uint64_t ith){
            info res;
            double first_left, first_right;
            if(m_s[ith] < m_t[ith]){
                res.split = source;
                if(ith == 0){
                    //Seed * (1+PathsFactorRight)
                    res.weight = m_s[ith] * (1+m_r[ith]);
                    res.first_left = false;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                first_right = m_s[ith] * ((1+m_r[ith]) + m_sol_r[ith] * (1+m_l[ith-1]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                first_left = m_s[ith] * ((1+m_l[ith-1]) + m_sol_l[ith-1] * (1+m_r[ith]));
            }else{
                res.split = target;
                if(ith == m_t.size()-1){
                    //Seed * (1+PathsFactorLeft)
                    res.weight = m_t[ith] * (1+m_l[ith]);
                    res.first_left = true;
                    return res;
                }
                //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
                first_right = m_t[ith] * ((1+m_r[ith+1]) + m_sol_r[ith+1] * (1+m_l[ith]));
                //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
                first_left = m_t[ith] * ((1+m_l[ith]) + m_sol_l[ith] * (1+m_r[ith+1]));
            }
            if(first_left <= first_right){
                res.weight = first_left;
                res.first_left = true;
            }else{
                res.weight = first_right;
                res.first_left = false;
            }
            return res;
        }

        info intersection(const uint64_t ith) {
            info res;
            res.split = intersect;
            double seed, first_right, first_left;
            seed = m_intersection[ith];
            //Seed * ((1+PathsFactorRight) + SolutionsFactorRight * (1+PathsFactorLeft))
            first_right = seed * ((1+m_r[ith+1]) + m_sol_r[ith+1] * (1+m_l[ith]));
            //Seed * ((1+PathsFactorLeft) + SolutionsFactorLeft * (1+PathsFactorRight))
            first_left = seed * ((1+m_l[ith]) + m_sol_l[ith] * (1+m_r[ith+1]));
            if(first_left <= first_right){
                res.weight = first_left;
                res.first_left = true;
            }else{
                res.weight = first_right;
                res.first_left = false;
            }
            return res;
        }

    };*/

    /*class h_distsigma_path {


    private:
        std::vector<double> m_s;
        std::vector<double> m_t;
        std::vector<double> m_r;
        std::vector<double> m_l;
        uint64_t m_max_p;
        uint64_t m_sigma;


    public:
        h_distsigma_path(const std::vector<PairPredPos> &preds,
                        const bwt_nose &L_S, const bwt_type &wt_pred_s,
                        uint64_t maxP, uint64_t sigma){

            m_max_p = maxP;
            m_sigma = sigma;
            //m_t.push_back(-1ULL);
            for(const auto &pair : preds){

                auto e_d = L_S.get_C(pair.id_pred + 1)-1;
                auto b_d = L_S.get_C(pair.id_pred);
                auto v_target = distinct_values(b_d, e_d, wt_pred_s);
                m_t.push_back(v_target / (double) m_sigma);

                auto rev_id = reverse(pair.id_pred, m_max_p);
                auto e_r = L_S.get_C(rev_id + 1)-1;
                auto b_r = L_S.get_C(rev_id);
                auto v_source = distinct_values(b_r, e_r, wt_pred_s);
                m_s.push_back(v_source / (double) m_sigma);
            }
            // m_s.push_back(-1ULL);
            auto s = m_s.size()+1;
            m_r.resize(s);
            m_r[s-1]=1;
            m_l.resize(s);
            m_l[0]=1;
            for(uint64_t i = 1; i < s; ++i){
                m_r[s-i-1] = m_r[s-i] * (1+m_t[m_t.size()-i]);
                m_l[i] = m_l[i-1] * (1+m_s[i-1]);
            }
            std::cout << "-----T-----" << std::endl;
            printVector(m_t);
            std::cout << "-----S-----" << std::endl;
            printVector(m_s);
            std::cout << "-----L-----" << std::endl;
            printVector(m_l);
            std::cout << "-----R-----" << std::endl;
            printVector(m_r);


        }

        info simple(const uint64_t ith){
            info res;
            //double a;
            if(m_s[ith] < m_t[ith]){
                res.split = source;
                if(ith == 0){
                    res.weight = m_s[ith] * m_r[ith+1];
                }else{
                    res.weight = m_s[ith] * m_r[ith+1] + m_s[ith] * m_l[ith-1];
                }
            }else{
                res.split = target;
                if(ith == m_t.size()-1){
                    res.weight = m_t[ith] * m_l[ith-1];
                }else{
                    res.weight = m_t[ith] * m_l[ith-1] + m_t[ith] * m_r[ith+1];
                }

            }
            return res;
        }

        info intersection(const uint64_t ith) {
            info res;
            res.split = intersect;
            double base, right, left;
            if(m_s[ith+1] < m_t[ith]){
                base = m_s[ith+1];
            }else{
                base = m_t[ith];
            }
            right = m_r[ith+2];
            left = m_l[ith];
            res.weight = left * base + base * right;
            return res;
        }

    };*/


}

#endif //RING_RPQ_SELECTIVITY_HPP
