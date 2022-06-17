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


namespace selectivity {

    enum split_type {target, source, intersect};

    inline uint64_t reverse(uint64_t id, uint64_t maxP){
        return (id > maxP) ? id - maxP : id + maxP;
    }

    inline uint64_t distinct_values(const uint64_t l, const uint64_t r, const bwt_type &wt_pred_s){
        auto p =  wt_pred_s.range_search_2d(l, r, 0, l-1, false);
        return p.first;
    }

    struct info {
        double  weight;
        split_type  split;
    };

    struct h_distinct {

        info simple(const uint64_t id,
                    const bwt_nose &L_S, const bwt_type &wt_pred_s,
                    uint64_t maxP, uint64_t sigma){
            auto e_d = L_S.get_C(id + 1)-1;
            auto b_d = L_S.get_C(id);
            auto rev_id = reverse(id, maxP);
            auto e_r = L_S.get_C(rev_id + 1)-1;
            auto b_r = L_S.get_C(rev_id);
            info res;
            auto v_target = distinct_values(b_d, e_d, wt_pred_s);
            auto v_source = distinct_values(b_r, e_r, wt_pred_s);
            std::cout << "Dist-source: " << v_source << std::endl;
            std::cout << "Dist-target: " << v_target << std::endl;
            if(v_source > v_target){
                res.split = target;
                res.weight = v_target / (double) sigma;
            }else{
                res.split = source;
                res.weight = v_source / (double) sigma;
            }
            return res;
        }

        info intersection(const uint64_t id1, const uint64_t id2,
                              const bwt_nose &L_S, const bwt_type &wt_pred_s,
                              uint64_t maxP, uint64_t sigma){
            auto e_l = L_S.get_C(id1 + 1)-1;
            auto b_l = L_S.get_C(id1);
            auto rev_id2 = reverse(id2, maxP);
            auto e_r = L_S.get_C(rev_id2 + 1)-1;
            auto b_r = L_S.get_C(rev_id2);
            double v_l = distinct_values(b_l, e_l, wt_pred_s);
            double v_r = distinct_values(b_r, e_r, wt_pred_s);
            std::cout << "Dist-left: " << v_l << std::endl;
            std::cout << "Dist-right: " << v_r << std::endl;
            info res;
            res.weight = (v_l * v_r) / (double) (sigma * sigma);
            res.split = intersect;
            return res;
        }
    };

    struct h_ratio {

        info simple(const uint64_t id,
                    const bwt_nose &L_S, const bwt_type &wt_pred_s,
                    uint64_t maxP, uint64_t sigma){
            auto e_d = L_S.get_C(id + 1)-1;
            auto b_d = L_S.get_C(id);
            auto rev_id = reverse(id, maxP);
            auto e_r = L_S.get_C(rev_id + 1)-1;
            auto b_r = L_S.get_C(rev_id);
            info res;
            auto v_source = distinct_values(b_d, e_d, wt_pred_s);
            auto v_target = distinct_values(b_r, e_r, wt_pred_s);
            if(v_source > v_target){
                res.split = target;
                res.weight = v_target / (double) (e_d-b_d+1);
            }else{
                res.split = source;
                res.weight = v_source / (double) (e_d-b_d+1);
            }
            return res;
        }

        info intersection(const uint64_t id1, const uint64_t id2,
                          const bwt_nose &L_S, const bwt_type &wt_pred_s,
                          uint64_t maxP, uint64_t sigma){
            auto e_l = L_S.get_C(id1 + 1)-1;
            auto b_l = L_S.get_C(id1);
            auto rev_id2 = reverse(id2, maxP);
            auto e_r = L_S.get_C(rev_id2 + 1)-1;
            auto b_r = L_S.get_C(rev_id2);
            auto v_l = distinct_values(b_l, e_l, wt_pred_s);
            auto v_r = distinct_values(b_r, e_r, wt_pred_s);
            info res;
            res.weight = (v_l / (double) (e_l-b_l+1)) * (v_r / (double) (e_r-b_r+1));
            res.split = intersect;
            return res;
        }
    };
}

#endif //RING_RPQ_SELECTIVITY_HPP
