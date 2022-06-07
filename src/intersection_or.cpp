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
#include <random>
#include <sdsl/int_vector.hpp>

uint64_t next(const std::vector<uint64_t>& vec, const uint64_t pos, const uint64_t value){
    for(uint64_t i = pos; i < vec.size(); ++i){
        if(vec[i]>=value){
            return i;
        }
    }
    return vec.size();

}

std::vector<uint64_t> base(8, 0);

typedef struct {
    uint64_t bytes;
    uint64_t codeword;
} cw_t;

cw_t encode(uint64_t i, uint64_t s){
    uint64_t c = 256-s;
    cw_t cw;
    cw.bytes = 1;
    cw.codeword = i % s;
    uint64_t x = i /s;
    while(x > 0){
        x = x -1;
        cw.codeword = (cw.codeword << 8);
        cw.codeword += ((x % c) + s);
        x = x / c;
        ++cw.bytes;
    }
    return cw;
}

 std::pair<uint64_t, uint64_t> decode(const sdsl::int_vector<> &buffer, uint64_t index, uint64_t s){
    uint64_t c = 256-s;
    uint64_t i = 0;
    uint64_t k = 0;
    while(buffer[index+k] >= s){
        i = i * c + (buffer[index+k]-s);
        ++k;
    }
    i = i * s + buffer[index+k];
    i = i + base[k];
    return {i, k};
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

   /* for(const auto &r : result){
       // std::cout << r << std::endl;
    }*/

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int32_t> dist(-20, 20);
    for(uint64_t i = 0; i < 20; ++i){
        std::cout << dist(gen) << std::endl;
    }

    uint64_t s = 230;
    base[0] = 0;
    for(uint64_t i = 1; i < 8; ++i){
        base[i] = base[i-1] + s * std::pow((256-s), (i-1));
    }



    std::vector<cw_t> cws;
    uint64_t bytes = 0;
    for(uint64_t i = 0; i < 100; ++i){
        auto cw = encode(i, s);
        cws.push_back(cw);
        bytes += cw.bytes;
        std::cout << "i: " << i << " cw.bytes:" << cw.bytes << " cw.codeword: " << cw.codeword << std::endl;
    }
    sdsl::int_vector<> buffer;
    buffer.width(8);
    buffer.resize(bytes);
    uint64_t index = 0;
    for(const auto &cw : cws){
        buffer.set_int(index, cw.codeword, cw.bytes*8);
        index += cw.bytes * 8;
    }
    index = 0;
    while(index < buffer.size()){
        auto pair = decode(buffer, index, s);
        std::cout << "Decoded: " << pair.first << std::endl;
        index = index + pair.second + 1;

    }





}