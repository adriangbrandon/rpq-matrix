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
// Created by Adrián on 26/5/22.
//

#ifndef RPQ_TREE_INCLUDED
#define RPQ_TREE_INCLUDED

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>
#include "regular.hpp"


/*typedef struct {
    int id_pred;
    int pos;

} PairPredPos;*/

typedef struct {
    int pred_b, pred_e;
    int pos_b, pos_e;

} PairPredPos;

typedef struct std::vector<PairPredPos> MandatoryData;

typedef struct {
    std::string pattern;
    std::unordered_map<char, uint64_t> id_to_pred;
    std::unordered_map<uint64_t, std::string> pred_to_str;
} PatternData;

class RpqTree {

private:

    void mandatoryPlusTraversal(Tree* e, MandatoryData &md, int& last);
    void mandatoryTraversal(Tree* e, MandatoryData &md, int& last);

    void splitTraversal(Tree* e, int split_pos, bool &left, std::string &rpq_l, std::string &rpq_r);

    std::string posToPredStr(int p);
    int posToPred(int pos);
public:

    RpqTree();

    RpqTree(const std::string &rpq, const std::unordered_map<std::string, uint64_t> &predicates, uint64_t max_pred);

    ~RpqTree(void);

    //TODO: debería devolver os predicados
    MandatoryData getMandatoryData();

    std::pair<std::string, std::string> splitRpq(int p);

    int patternPredicates();

    Tree* root();

    int getPred(int pos);




private:
    Tree *tree;
    Tree **pos;
    int m;
    int* pos_id;
    PatternData patternData;
};

#endif //RPQ_TREE_INCLUDED
