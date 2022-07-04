/*
 * bwt_interval.hpp
 * Copyright (C) 2020 Author removed for double-blind evaluation
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

#ifndef BWT_INTERVAL
#define BWT_INTERVAL

#include "bwt.hpp"

class bwt_interval 
{
    uint64_t l;
    uint64_t r;
    uint64_t cur_val;  // current value within the interval
    uint64_t cur_rank;
public: 
    bwt_interval() {;}


    ~bwt_interval() {;}


    bwt_interval(uint64_t _l, uint64_t _r)
    {
        l = _l;
        r = _r;
        cur_val = (uint64_t)-1;
        cur_rank = (uint64_t)-1;
    }


    bwt_interval(const bwt_interval & i)
    {
        l = i.l;
        r = i.r;
        cur_val = i.cur_val;
        cur_rank = i.cur_rank;
    }

    uint64_t begin(bwt & B)
    {
        return B.min_in_range(l, r);
    }

    uint64_t next_value(uint64_t val, bwt & B)
    {
        return B.range_next_value(val, l, r);
    }

    inline uint64_t end()
    {
        return 0;  // assumes 0 is not stored in the triples
    }

    inline uint64_t get_cur_value()
    {
        return cur_val;
    }

    inline std::pair<uint64_t, uint64_t>
    get_stored_values()
    {
        return std::pair<uint64_t, uint64_t>(cur_val, cur_rank);
    }

    inline void 
    set_stored_values(uint64_t _cur_val, uint64_t _r)
    {
        cur_val = _cur_val;
        cur_rank = _r;
    }

    inline void set_right(uint64_t _r){
        r = _r;
    }

    bwt_interval operator =(bwt_interval i)
    {
        l = i.l;
        r = i.r;
        cur_val = i.cur_val;
        cur_rank = i.cur_rank;
        return *this;
    }

    inline uint64_t left() const
    {
        return l;
    }

    inline uint64_t right() const
    {
        return r;
    }

    inline uint64_t size()
    {
        return r-l+1;
    }
};

#endif
