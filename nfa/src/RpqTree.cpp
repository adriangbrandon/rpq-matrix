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

#include "RpqTree.hpp"


PatternData rpqToPatternData(const std::string &rpq,
                             const std::unordered_map<std::string, uint64_t> &predicates,
                             uint64_t max_pred)
{
    //std::unordered_map<uint64_t, char> map_pred_to_id;
    std::unordered_map<char, uint64_t> map_id_to_pred;
    std::unordered_map<uint64_t , std::string> map_pred_to_str;
    char firstUnassigned = 1;
    auto pattern = rpq;

    std::regex predRE("(<[^>]*>)");
    std::sregex_iterator next(rpq.begin(), rpq.end(), predRE);
    std::sregex_iterator end;

    while (next != end)
    {
        std::smatch match = *next;
        const auto str = match.str();
        uint64_t predicateId;
        if(str[1] == '%'){ //<%...>
            std::string no_neg = '<' + str.substr(2);
            predicateId = predicates.at(no_neg) + max_pred;
        }else{
            predicateId = predicates.at(str);
        }
        if (map_pred_to_str.find(predicateId) == map_pred_to_str.end())
        {
            map_pred_to_str[predicateId] = str;
            map_id_to_pred[firstUnassigned] = predicateId;

            std::regex replaceRE(str);
            pattern = std::regex_replace(pattern, replaceRE, std::string(1, firstUnassigned));

            ++firstUnassigned;
        }
        next++;
    }

    std::regex slashRE("/");
    pattern = std::regex_replace(pattern, slashRE, "");

    return {pattern, map_id_to_pred, map_pred_to_str};
}

void RpqTree::mandatoryPlusTraversal(Tree *e, MandatoryData &md, int &last) {

    Tree* left_node = e->e1;
    Tree* right_node = e->e2;
    while(left_node->type == CONC || left_node->type == PLUS){
        left_node = left_node->e1;
    }
    while(right_node->type == CONC || right_node->type == PLUS){
        if(right_node->type == CONC){
            right_node = right_node->e2;
        }else{
            right_node = right_node->e1;
        }
    }
    if(left_node->type == STR && right_node->type == STR) {
        md.push_back({posToPred(left_node->pos), posToPred(right_node->pos), left_node->pos, right_node->pos});
        last = right_node->pos;
    }
}

void RpqTree::mandatoryTraversal(Tree* e, MandatoryData &md, int& last){
    if(e->type == STR){
        auto pred = posToPred(e->pos);
        md.push_back({pred, pred, e->pos, e->pos});
        last = e->pos;
    }else if(e->type == CONC){
        mandatoryTraversal(e->e1, md, last);
        mandatoryTraversal(e->e2, md, last);
    }else if(e->type == PLUS){
        if(e->e1->type == STR){
            auto pred = posToPred(e->e1->pos);
            md.push_back({pred, pred, e->e1->pos, e->e1->pos});
            last = e->e1->pos;
        }else if (e->e1->type == CONC){
            mandatoryPlusTraversal(e->e1,  md, last);
        }

    }
}

void RpqTree::splitTraversal(Tree* e, int split_pos, bool &left, std::string &rpq_l, std::string &rpq_r){

    switch (e->type)
    {
        case STR:
            if(left){
                rpq_l += posToPredStr(e->pos);
                if(e->pos == split_pos) left = false;
            }else{
                rpq_r += posToPredStr(e->pos);
            }
            break;
        case STAR:
            if(left) {
                rpq_l += "(";
                splitTraversal(e->e1, split_pos, left, rpq_l, rpq_r);
                rpq_l += ")*";
            }else{
                rpq_r += "(";
                splitTraversal(e->e1, split_pos, left, rpq_l, rpq_r);
                rpq_r += ")*";
            }
            break;
        case PLUS:
            if(left) {
                rpq_l += "(";
                splitTraversal(e->e1, split_pos, left, rpq_l, rpq_r);
                rpq_l += ")+";
            }else{
                rpq_r += "(";
                splitTraversal(e->e1, split_pos, left, rpq_l, rpq_r);
                rpq_r += ")+";
            }
            break;
        case QUESTION:
            if(left) {
                rpq_l += "(";
                splitTraversal(e->e1, split_pos, left, rpq_l, rpq_r);
                rpq_l += ")?";
            }else{
                rpq_r += "(";
                splitTraversal(e->e1, split_pos, left, rpq_l, rpq_r);
                rpq_r += ")?";
            }
            break;
        case OOR:
            if(left){
                rpq_l += "(";
                splitTraversal(e->e1, split_pos, left, rpq_l, rpq_r);
                rpq_l += "|";
                splitTraversal(e->e2, split_pos, left, rpq_l, rpq_r);
                rpq_l += ")";
            }else{
                rpq_r += "(";
                splitTraversal(e->e1, split_pos, left, rpq_l, rpq_r);
                rpq_r += "|";
                splitTraversal(e->e2, split_pos, left, rpq_l, rpq_r);
                rpq_r += ")";
            }
            break;
        case CONC:
            splitTraversal(e->e1, split_pos, left, rpq_l, rpq_r);
            splitTraversal(e->e2, split_pos, left, rpq_l, rpq_r);
            break;

    }
}

inline int RpqTree::posToPred(int p){
    return patternData.id_to_pred[pos_id[p]];
}

inline std::string RpqTree::posToPredStr(int p){
    return patternData.pred_to_str[patternData.id_to_pred[pos_id[p]]];
}

RpqTree::RpqTree(const std::string &rpq, const std::unordered_map<std::string, uint64_t> &predicates,
                 uint64_t max_pred){

    patternData = rpqToPatternData(rpq, predicates, max_pred);
    m = patternData.pattern.length();
    pos = new Tree *[m];
    pos_id = new int[m+1];
    tree = parse(patternData.pattern.c_str(), m, pos, false);
    regularPreprocBInv(patternData.pattern.c_str(), tree, pos, pos_id);
}

MandatoryData RpqTree::getMandatoryData() {
    MandatoryData md;
    int last = 0;
    mandatoryTraversal(tree, md, last);
    return md;
}

std::pair<std::string, std::string> RpqTree::splitRpq( int p_split){
    std::pair<std::string, std::string> res;
    bool left = (p_split > 0);
    splitTraversal(tree, p_split, left, res.first,  res.second);
    return res;
}

int RpqTree::patternPredicates() {
    return patternData.id_to_pred.size();
}

int RpqTree::getPred(int p){
    return patternData.id_to_pred[pos_id[p]];
}

Tree *RpqTree::root() {
    return this->tree;
}


/*std::pair<std::string, std::string> RpqTree::splitRpq(const std::string &rpq, int p_split){
    int i = 0, j = 0, p = 0;
    while(i < m){
        if(pos[i] != NULL){
            ++i;
            j = j + pos_length[p];
            ++p;
            if(pos[i-1]->pos == p_split){
                while(i < m && (pos[i] == NULL && patternData.pattern[i] != '(' && patternData.pattern[i] != '|')){
                    ++i;
                    ++j;
                }
                break;
            }
        }else{
            ++i;
            ++j;
        }
    }
    std::pair<std::string, std::string> res;
    res.first = rpq.substr(0, j);
    res.second = rpq.substr(j);
    return res;
}*/



RpqTree::~RpqTree() {
    freeTree(tree);
    delete [] pos;
    delete [] pos_id;
}
