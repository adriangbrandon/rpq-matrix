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
// Created by Adrián on 14/5/22.
//

#include <string>
#include <stack>
#include <iostream>
#include <vector>

int main(int argc, char **argv) {

    if(argc != 3){
        std::cout << argv[0] << " <rpq> <n_preds>" << std::endl;
    }
    std::string rpq = argv[1];
    uint64_t n_preds = atoi(argv[2]);

    typedef struct {
        uint64_t beg;
        bool is_or = false;
    } section_type;

    typedef struct {
        uint64_t beg;
        uint64_t end;
        bool mandatory = true;
    } pred_data_type;

    uint64_t start;
    uint64_t i = 0, p_i;
    std::vector<pred_data_type> pred_v;
    std::stack<section_type> sections; //first + parenthesis
    section_type first;
    first.beg = 0;
    sections.push(first);
    for(uint64_t p = 0; p < n_preds; ++p){
        while(rpq.at(i) == '(') {
            section_type sec;
            sec.beg = p;
            sections.push(sec);
            ++i;
        }
        pred_data_type pred_data;
        pred_data.beg = i;
        //3. Find the end of predicate
        while (rpq.at(i) != '>') ++i;
        ++i; //Next to '>'
        pred_data.end = i;
        pred_v.push_back(pred_data);
        while(i < rpq.size() && (rpq.at(i)!='<' && rpq.at(i)!='(')){
            if (rpq.at(i)==')') {
                auto sec = sections.top();
                sections.pop();
                p_i = sec.beg;
                if (sec.is_or) {
                    for(auto j = p_i; j <= p; ++j) {
                        pred_v[j].mandatory = false;
                    }
                }
            }else if(rpq.at(i)=='|'){
                sections.top().is_or=true;
            }else if (rpq.at(i)=='?' || rpq.at(i)=='*'){
                for(auto j = p_i; j <= p; ++j) { //p_i is initialized because it has to be preceded by ')'
                    pred_v[j].mandatory = false;
                }
            }
            ++i;
        }
        if(i == rpq.size()){
            if(sections.top().is_or) { //first section (no parentheses)
                auto sec = sections.top();
                sections.pop();
                p_i = sec.beg;
                if (sec.is_or) {
                    for(auto j = p_i; j <= p; ++j) {
                        pred_v[j].mandatory = false;
                    }
                }
            }
        }

    }

    for(const auto &p : pred_v){
        std::cout << rpq.substr(p.beg, p.end - p.beg) << " mandatory: " << (uint64_t ) p.mandatory << std::endl;
    }



}
