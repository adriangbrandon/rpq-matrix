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
// Created by Adrián on 17/5/22.
//

#ifndef RING_RPQ_WT_INTERSECTION_HPP
#define RING_RPQ_WT_INTERSECTION_HPP

#include <algorithm>
#include <utility>
#include <sdsl/wt_helper.hpp>

namespace sdsl {

    template<class t_wt>
    std::vector<typename t_wt::value_type>
    intersect_nofreq(const t_wt& wt, const std::vector<range_type>& ranges, typename t_wt::size_type t=0)
    {
        using std::get;
        using size_type      = typename t_wt::size_type;
        using value_type     = typename t_wt::value_type;
        using node_type      = typename t_wt::node_type;
        using pnvr_type      = std::pair<node_type, range_vec_type>;
        typedef std::stack<pnvr_type> stack_type;

        static_assert(has_expand<t_wt, std::array<node_type,2>(const node_type&)>::value,
                "intersect requires t_wt to have expand(const node_type&)");

        std::vector<size_type> res;

        auto push_node = [&wt,&t](stack_type& s, node_type& child,
                                  range_vec_type& child_range) {
            auto end = std::remove_if(child_range.begin(), child_range.end(),
                                      [&](const range_type& x) { return empty(x);});
            if (end > child_range.begin() + t - 1) {
                s.emplace(pnvr_type(child, range_vec_type(child_range.begin(),
                                                          end)));
            }
        };

        if (ranges.empty())
            return res;

        t = (t==0) ? ranges.size() : t;

        std::stack<pnvr_type> stack;
        stack.emplace(pnvr_type(wt.root(), ranges));

        while (!stack.empty()) {
            pnvr_type x = stack.top(); stack.pop();

            if (wt.is_leaf(x.first)) {
                const auto& iv = x.second;
                if (t <= iv.size()) {
                    res.emplace_back(wt.sym(x.first));
                }
            } else {
                auto child        = wt.expand(x.first);
                auto child_ranges = wt.expand(x.first, x.second);

                push_node(stack, get<1>(child), get<1>(child_ranges));
                push_node(stack, get<0>(child), get<0>(child_ranges));
            }
        }
        return res;
    }

}

#endif //RING_RPQ_WT_INTERSECTION_HPP
