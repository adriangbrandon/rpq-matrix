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
// Created by Adrián on 16/5/22.
//


#include <vector>
#include <iostream>

uint64_t next(const std::vector<uint64_t>& vec, const uint64_t pos, const uint64_t value){
    for(uint64_t i = pos; i < vec.size(); ++i){
        if(vec[i]>=value){
            return i;
        }
    }
    return vec.size();

}


int main(int argc, char **argv) {

    std::vector<uint64_t> v1 = {10, 20, 30, 40, 50, 60, 70};
    std::vector<uint64_t> or1 = {2, 4, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
    std::vector<uint64_t> or2 = {30, 35, 40, 45, 55, 65, 75};
    std::vector<uint64_t> or3 = {70,80,90};
    std::vector<std::vector<uint64_t>> or_all = {or1, or2, or3};
    std::vector<uint64_t> or_pos(or_all.size(), 0);

    std::vector<uint64_t> result;
    uint64_t value = 0, pos1 = 0;
    while(pos1 < v1.size()){
        pos1 = next(v1, pos1, value);
        value = v1[pos1];
        for(uint64_t j = 0; j < or_all.size(); ++j){
            or_pos[j] = next(or_all[j], or_pos[j], value);
        }
        bool ok = false;
        uint64_t min = -1ULL;
        for(uint64_t j = 0; j < or_all.size(); ++j){
            if(or_pos[j] == or_all[j].size()) continue;
            uint64_t v_or = or_all[j][or_pos[j]];
            if(value == v_or) ok = true;
            if(v_or < min) min = v_or;
        }
        if(ok) {
            result.push_back(value);
            ++value;
        }else{
            value = min;
        }
    }

    for(const auto &r : result){
        std::cout << r << std::endl;
    }

}