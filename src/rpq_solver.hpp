//
// Created by Adrián on 3/5/23.
//

#ifndef MATRIX_RPQ_RPQ_SOLVER_HPP
#define MATRIX_RPQ_RPQ_SOLVER_HPP

#define N 958844164
#define SIZE 5420 // 1 to 5419
#define V 296008192 // 1 to...

extern "C" {
#include <matrix.h>
}


#include <RpqTree.hpp>
#include <list>
#include <vector>

namespace rpq {

    class solver {

        typedef struct {
            matrix m;
            bool is_tmp;
            bool is_fixed;
        } data_type;
        typedef std::list<data_type> list_type;
        typedef typename list_type::iterator it_type;
        typedef typename list_type::reverse_iterator rev_it_type;

    private:

        std::vector<matrix> m_matrices;
        //std::vector<matrix> m_tmp_matrices;
        std::unordered_map<std::string, uint64_t> m_map_SO;
        std::unordered_map<std::string, uint64_t> m_map_P;

        void traversal(RpqTree* rpqTree, Tree* node, int parentType, list_type &res,
                       bool skip_closure = false){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    res.insert(res.begin(), data_type{m_matrices[pred], false, false});
                    break;
                }
                case CONC:
                {
                    list_type ll, rl;
                    traversal(rpqTree, node->e1, CONC, ll);
                    traversal(rpqTree, node->e2, CONC, rl);
                    rl.splice(rl.begin(), ll);
                    if (parentType == CONC){
                        res = std::move(rl);
                    }else{
                        it_type it1, it2, it1_min, it2_min;
                        matrix tmp;
                        while(rl.size() > 2){//Check size if there are more than three elements
                            it1 = it2 = rl.begin();
                            ++it2;
                            uint64_t min = UINT64_MAX, weight;
                            while(it2 != rl.end()){
                                weight = matSpace(it1->m) + matSpace(it2->m);
                                if(min > weight){
                                    it1_min = it1;
                                    it2_min = it2;
                                    min = weight;
                                }
                                ++it1; ++it2;
                            }
                            tmp = matMult(it1_min->m, it2_min->m);
                            if(it1_min->is_tmp) matDestroy(it1_min->m);
                            if(it2_min->is_tmp) matDestroy(it2_min->m);

                            rl.insert(it1_min, data_type{tmp, true, false});
                            rl.erase(it1_min);
                            rl.erase(it2_min);
                        }
                        it1_min = it2_min = rl.begin();
                        ++it2_min;
                        tmp = matMult(it1_min->m, it2_min->m);
                        if(it1_min->is_tmp) matDestroy(it1_min->m);
                        if(it2_min->is_tmp) matDestroy(it2_min->m);
                        res.insert(res.begin(), data_type{tmp, true, false});
                    }
                    break;
                }

                case OOR:
                {
                    list_type ll, rl;
                    traversal(rpqTree, node->e1, OOR, ll);
                    traversal(rpqTree, node->e2, OOR, rl);
                    rl.splice(rl.begin(), ll);
                    if (parentType == OOR){
                        res = std::move(rl);
                    }else{
                        it_type it1, it2, it1_min, it2_min;
                        matrix tmp;
                        while(rl.size() > 2){//Check size if there are more than three elements
                            it1 = it2 = rl.begin();
                            uint64_t min = UINT64_MAX, weight;
                            while(it1 != rl.end()){
                                weight = matSpace(it1->m);
                                if(min > weight){
                                    it1_min = it1;
                                    min = weight;
                                }
                                ++it1;
                            }
                            min = UINT64_MAX;
                            while(it2 != rl.end()){
                                if(it2 != it1_min ){
                                    weight = matSpace(it2->m);
                                    if(min > weight){
                                        it2_min = it2;
                                        min = weight;
                                    }
                                }
                                ++it2;
                            }
                            tmp = matSum(it1_min->m, it2_min->m);
                            if(it1_min->is_tmp) matDestroy(it1_min->m);
                            if(it2_min->is_tmp) matDestroy(it2_min->m);

                            rl.insert(it1_min, data_type{tmp, true, false});
                            rl.erase(it1_min);
                            rl.erase(it2_min);
                        }

                        it1_min = it2_min = rl.begin();
                        ++it2_min;
                        tmp = matSum(it1_min->m, it2_min->m);
                        if(it1_min->is_tmp) matDestroy(it1_min->m);
                        if(it2_min->is_tmp) matDestroy(it2_min->m);
                        res.insert(res.begin(), data_type{tmp, true, false});
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, STAR, ll);
                    if(!skip_closure){
                        matrix tmp = matClos(ll.front().m, 0);
                        if(ll.front().is_tmp) matDestroy(ll.front().m);
                        // std::cout << "STAR : " << ll.front() << std::endl;
                        res.insert(res.begin(), data_type{tmp, true, false});
                    }else{
                        res = std::move(ll);
                    }
                    break;
                }
                case PLUS:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, PLUS, ll);
                    if(!skip_closure){
                        matrix tmp = matClos(ll.front().m, 1);
                        if(ll.front().is_tmp) matDestroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, true, false});
                    }else{
                        res = std::move(ll);
                    }
                    break;
                }
                case QUESTION:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, QUESTION, ll);
                    matrix Id = matId(std::max(ll.front().m->height,ll.front().m->width));
                    matrix tmp = matSum(ll.front().m, Id);
                    matDestroy(Id);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true, false});
                    break;
                }
            }
        }

        void traversal_col_fixed(RpqTree* rpqTree, Tree* node, int parentType, int col, list_type &res,
                                 bool skip_closure = false){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    matrix a = m_matrices[pred];
                    if(parentType != ROOT){
                        res.insert(res.begin(), data_type{a, false, false});
                    }else {
                        matrix e = matEmpty(a->height, a->width);
                        matrix m = matSum1(fullSide, a, e, col);
                        matDestroy(e);
                        res.insert(res.begin(), data_type{m, true, true});
                    }
                    break;
                }
                case CONC:
                {
                    list_type ll, rl;
                    if(node->e1->type == STAR || node->e1->type == PLUS){ //Special case!
                        traversal(rpqTree, node->e1, CONC, ll, true);
                        traversal_col_fixed(rpqTree, node->e2, CONC, col, rl);
                        matrix tmp;
                        uint pos = (node->e1->type == STAR) ? 0 : 1;
                        if(rl.size() > 1){ //product of the right part
                            auto it2 = rl.rbegin(); //last element
                            auto it1 = rl.rbegin(); //last element
                            ++it2; //last element -1
                            matrix aux;
                            if(it1->is_fixed){
                                tmp = matMult(it2->m, it1->m);
                            }else{
                                tmp = matMult1(fullSide, it2->m,it1->m, col);
                            }
                            if(it1->is_tmp) matDestroy(it1->m);
                            if(it2->is_tmp) matDestroy(it2->m);
                            while(++it2 != rl.rend()){
                                aux = matMult(it2->m, tmp); //tmp is already fixed
                                if(it2->is_tmp) matDestroy(it2->m);
                                matDestroy(tmp);
                                tmp = aux;
                            }
                            aux = matClosMult1(fullSide, ll.front().m, pos, tmp, col);
                            matDestroy(tmp);
                            if(ll.front().is_tmp) matDestroy(ll.front().m);
                            res.insert(res.begin(), data_type{aux, true, true});
                        }else{
                            tmp = matClosMult1(fullSide, ll.front().m, pos, rl.front().m, col);
                            if(ll.front().is_tmp) matDestroy(ll.front().m);
                            if(rl.front().is_tmp) matDestroy(rl.front().m);
                            res.insert(res.begin(), data_type{tmp, true, true});
                        }
                    }else{
                        traversal(rpqTree, node->e1, CONC, ll);
                        traversal_col_fixed(rpqTree, node->e2, CONC, col, rl);
                        rl.splice(rl.begin(), ll);
                        if (parentType == CONC){
                            res = std::move(rl);
                        }else{
                            auto it2 = rl.rbegin(); //last element
                            auto it1 = rl.rbegin(); //last element
                            ++it2; //last element -1
                            matrix tmp, aux;
                            if(it1->is_fixed){
                                tmp = matMult(it2->m, it1->m);
                            }else{
                                tmp = matMult1(fullSide, it2->m,it1->m, col);
                            }
                            if(it1->is_tmp) matDestroy(it1->m);
                            if(it2->is_tmp) matDestroy(it2->m);
                            while(++it2 != rl.rend()){
                                aux = matMult(it2->m, tmp); //tmp is already fixed
                                if(it2->is_tmp) matDestroy(it2->m);
                                matDestroy(tmp);
                                tmp = aux;
                            }
                            res.insert(res.begin(), data_type{tmp, true, true});
                        }
                    }
                    break;
                }

                case OOR:
                {
                    list_type ll, rl;
                    traversal_col_fixed(rpqTree, node->e1, OOR, col, ll);
                    traversal_col_fixed(rpqTree, node->e2, OOR, col, rl);
                    rl.splice(rl.begin(), ll);
                    if (parentType == OOR){
                        res = std::move(rl);
                    }else{
                        matrix tmp;
                        it_type it1, it2, it1_min, it2_min;
                        while(rl.size() > 2) {//Check if there a
                            it1 = it2 = rl.begin();
                            uint64_t min = UINT64_MAX, weight;
                            while (it1 != rl.end()) {
                                weight = matSpace(it1->m);
                                if (min > weight) {
                                    it1_min = it1;
                                    min = weight;
                                }
                                ++it1; //++it2;
                            }
                            min = UINT64_MAX;
                            while (it2 != rl.end()) {
                                if (it2 != it1_min) {
                                    weight = matSpace(it2->m);
                                    if (min > weight) {
                                        it2_min = it2;
                                        min = weight;
                                    }
                                }
                                ++it2;
                            }
                            //tmp = matSum1(it1_min->m, it2_min->m, fullSide, col);
                            if(it1_min->is_fixed && it2_min->is_fixed){
                                tmp = matSum(it1_min->m, it2_min->m);
                            }else{
                                tmp = matSum1(fullSide, it1_min->m, it2_min->m, col);
                            }
                            if (it1_min->is_tmp) matDestroy(it1_min->m);
                            if (it2_min->is_tmp) matDestroy(it2_min->m);

                            rl.insert(it1_min, data_type{tmp, true, true});
                            rl.erase(it1_min);
                            rl.erase(it2_min);
                        }
                        it1_min = it2_min = rl.begin();
                        ++it2_min;
                        if(it1_min->is_fixed && it2_min->is_fixed){
                            tmp = matSum(it1_min->m, it2_min->m);
                        }else{
                            tmp = matSum1(fullSide, it1_min->m, it2_min->m, col);
                        }
                        if(it1_min->is_tmp) matDestroy(it1_min->m);
                        if(it2_min->is_tmp) matDestroy(it2_min->m);
                        res.insert(res.begin(), data_type{tmp, true, true});
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, STAR, ll);
                    if(!skip_closure){
                        matrix tmp = matClos1(fullSide, ll.front().m, 0, col);
                        // std::cout << "STAR : " << ll.front() << std::endl;
                        if(ll.front().is_tmp) matDestroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, true, true});
                    }else{
                        res = std::move(ll);
                    }

                    break;
                }
                case PLUS:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, PLUS, ll);
                    if(!skip_closure){
                        matrix tmp = matClos1(fullSide, ll.front().m, 1, col);
                        if(ll.front().is_tmp) matDestroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, true, true});
                    }else{
                        res = std::move(ll);
                    }
                    break;
                }
                case QUESTION:
                {
                    list_type ll;
                    traversal_col_fixed(rpqTree, node->e1, QUESTION, col, ll);
                    matrix Id = matId(std::max(ll.front().m->height,ll.front().m->width));
                    matrix tmp = matSum1(fullSide, ll.front().m, Id, col);
                    matDestroy(Id);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true, true});
                    break;
                }
            }
        }

        void traversal_row_fixed(RpqTree* rpqTree, Tree* node, int parentType, int row, list_type &res,
                                 bool skip_closure = false){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    matrix a = m_matrices[pred];
                    if(parentType != ROOT){
                        res.insert(res.begin(), data_type{a, false, false});
                    }else {
                        matrix e = matEmpty(a->height, a->width);
                        matrix m = matSum1(row, a, e, fullSide);
                        matDestroy(e);
                        res.insert(res.begin(), data_type{m, true, true});
                    }
                    break;
                }
                case CONC:
                {
                    list_type ll, rl;
                    if(node->e2->type == STAR || node->e2->type == PLUS) { //Special case!
                        traversal_row_fixed(rpqTree, node->e1, CONC, row,ll);
                        traversal(rpqTree, node->e2, CONC, rl, true);
                        matrix tmp;
                        uint pos = (node->e2->type == STAR) ? 0 : 1;
                        if (ll.size() > 1) { //product of the right part
                            auto it2 = ll.begin(); //first element
                            auto it1 = ll.begin(); //first element
                            ++it2; //first element +1
                            matrix aux;
                            if(it1->is_fixed){
                                tmp = matMult(it1->m, it2->m);
                            }else{
                                tmp = matMult1(row, it1->m,it2->m, fullSide);
                            }
                            if(it1->is_tmp) matDestroy(it1->m);
                            if(it2->is_tmp) matDestroy(it2->m);
                            while(++it2 != ll.end()){
                                aux = matMult(it2->m, tmp); //tmp is already fixed
                                if(it2->is_tmp) matDestroy(it2->m);
                                matDestroy(tmp);
                                tmp = aux;
                            }
                            aux = matMultClos1(row, tmp, rl.front().m, pos, fullSide);
                            matDestroy(tmp);
                            if(rl.front().is_tmp) matDestroy(rl.front().m);
                            res.insert(res.begin(), data_type{aux, true, true});
                        }else{
                            tmp = matMultClos1(row, ll.front().m, rl.front().m, pos, fullSide);
                            if(ll.front().is_tmp) matDestroy(ll.front().m);
                            if(rl.front().is_tmp) matDestroy(rl.front().m);
                            res.insert(res.begin(), data_type{tmp, true, true});
                        }
                    }else{
                        traversal_row_fixed(rpqTree, node->e1, CONC, row,ll);
                        traversal(rpqTree, node->e2, CONC, rl);
                        rl.splice(rl.begin(), ll);
                        if (parentType == CONC){
                            res = std::move(rl);
                        }else{
                            auto it2 = rl.begin(); //first element
                            auto it1 = rl.begin(); //first element
                            ++it2; //first element +1
                            matrix tmp, aux;
                            if(it1->is_fixed){
                                tmp = matMult(it1->m, it2->m);
                            }else{
                                tmp = matMult1(row, it1->m,it2->m, fullSide);
                            }
                            if(it1->is_tmp) matDestroy(it1->m);
                            if(it2->is_tmp) matDestroy(it2->m);
                            while(++it2 != rl.end()){
                                aux = matMult(tmp, it2->m);
                                if(it2->is_tmp) matDestroy(it2->m);
                                matDestroy(tmp);
                                tmp = aux;
                            }
                            res.insert(res.begin(), data_type{tmp, true, true});
                        }
                    }
                    break;
                }
                case OOR:
                {
                    list_type ll, rl;
                    traversal_row_fixed(rpqTree, node->e1, OOR, row, ll);
                    traversal_row_fixed(rpqTree, node->e2, OOR, row, rl);
                    rl.splice(rl.begin(), ll);
                    if (parentType == OOR){
                        res = std::move(rl);
                    }else{

                        matrix tmp;
                        it_type it1, it2, it1_min, it2_min;
                        while(rl.size() > 2) {//Check if there a
                            it1 = it2 = rl.begin();
                            uint64_t min = UINT64_MAX, weight;
                            while (it1 != rl.end()) {
                                weight = matSpace(it1->m);
                                if (min > weight) {
                                    it1_min = it1;
                                    min = weight;
                                }
                                ++it1; //++it2;
                            }
                            min = UINT64_MAX;
                            while (it2 != rl.end()) {
                                if (it2 != it1_min) {
                                    weight = matSpace(it2->m);
                                    if (min > weight) {
                                        it2_min = it2;
                                        min = weight;
                                    }
                                }
                                ++it2;
                            }
                            //tmp = matSum1(it1_min->m, it2_min->m, row, fullSide);
                            if(it1_min->is_fixed && it2_min->is_fixed){
                                tmp = matSum(it1_min->m, it2_min->m);
                            }else{
                                tmp = matSum1(row, it1_min->m, it2_min->m, fullSide);
                            }
                            if (it1_min->is_tmp) matDestroy(it1_min->m);
                            if (it2_min->is_tmp) matDestroy(it2_min->m);

                            rl.insert(it1_min, data_type{tmp, true, true});
                            rl.erase(it1_min);
                            rl.erase(it2_min);
                        }
                        it1_min = it2_min = rl.begin();
                        ++it2_min;
                        //tmp = matSum1(it1_min->m, it2_min->m, row, fullSide);
                        if(it1_min->is_fixed && it2_min->is_fixed){
                            tmp = matSum(it1_min->m, it2_min->m);
                        }else{
                            tmp = matSum1(row, it1_min->m, it2_min->m, fullSide);
                        }
                        if(it1_min->is_tmp) matDestroy(it1_min->m);
                        if(it2_min->is_tmp) matDestroy(it2_min->m);
                        res.insert(res.begin(), data_type{tmp, true, true});
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, STAR, ll);
                    if(!skip_closure){
                        matrix tmp = matClos1(row, ll.front().m, 0, fullSide);
                        // std::cout << "STAR : " << ll.front() << std::endl;
                        if(ll.front().is_tmp) matDestroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, true, true});
                    }else{
                        res = std::move(ll);
                    }
                    break;
                }
                case PLUS:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, PLUS, ll);
                    if(!skip_closure){
                        matrix tmp = matClos1(row, ll.front().m, 1, fullSide);
                        if(ll.front().is_tmp) matDestroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, true, true});
                        break;
                    }else{
                        res = std::move(ll);
                    }
                }
                case QUESTION:
                {
                    list_type ll;
                    traversal_row_fixed(rpqTree, node->e1, QUESTION, row, ll);
                    matrix Id = matId(std::max(ll.front().m->height,ll.front().m->width));
                    matrix tmp = matSum1(row, ll.front().m, Id, fullSide);
                    matDestroy(Id);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true, true});
                    break;
                }
            }
        }

        void traversal_row_col_fixed(RpqTree* rpqTree, Tree* node, int parentType, int row, int col, list_type &res,
                                     bool skip_closure = false){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    matrix a = m_matrices[pred];
                    if(parentType != ROOT){
                        res.insert(res.begin(), data_type{a, false, false});
                    }else {
                        matrix e = matEmpty(a->height, a->width);
                        matrix m = matSum1(row, a, e, col);
                        res.insert(res.begin(), data_type{m, true, true});
                        matDestroy(e);
                    }
                    break;
                }
                case CONC:
                {
                    list_type ll, rl;
                    traversal_row_fixed(rpqTree, node->e1, CONC, row,ll);
                    traversal_col_fixed(rpqTree, node->e2, CONC, col, rl);
                    rl.splice(rl.begin(), ll);
                    if (parentType == CONC){
                        res = std::move(rl);
                    }else{
                        auto it2 = rl.begin(); //first element
                        auto it1 = rl.begin(); //first element
                        ++it2; //first element +1
                        matrix tmp, aux;
                        if(it1->is_fixed){
                            tmp = matMult(it1->m, it2->m);
                        }else{
                            tmp = matMult1(row, it1->m,it2->m, col);
                        }
                        if(it1->is_tmp) matDestroy(it1->m);
                        if(it2->is_tmp) matDestroy(it2->m);
                        while(++it2 != rl.end()){
                            aux = matMult(it2->m,tmp);
                            if(it2->is_tmp) matDestroy(it2->m);
                            matDestroy(tmp);
                            tmp = aux;
                        }
                        res.insert(res.begin(), data_type{tmp, true, true});
                    }
                    break;
                }
                case OOR:
                {
                    list_type ll, rl;
                    traversal_row_col_fixed(rpqTree, node->e1, OOR, row, col, ll);
                    traversal_row_col_fixed(rpqTree, node->e2, OOR, row, col, rl);
                    rl.splice(rl.begin(), ll);
                    if (parentType == OOR){
                        res = std::move(rl);
                    }else{

                        matrix tmp;
                        it_type it1, it2, it1_min, it2_min;
                        while(rl.size() > 2) {//Check if there a
                            it1 = it2 = rl.begin();
                            uint64_t min = UINT64_MAX, weight;
                            while (it1 != rl.end()) {
                                weight = matSpace(it1->m);
                                if (min > weight) {
                                    it1_min = it1;
                                    min = weight;
                                }
                                ++it1; //++it2;
                            }
                            min = UINT64_MAX;
                            while (it2 != rl.end()) {
                                if (it2 != it1_min) {
                                    weight = matSpace(it2->m);
                                    if (min > weight) {
                                        it2_min = it2;
                                        min = weight;
                                    }
                                }
                                ++it2;
                            }
                            //tmp = matSum1(it1_min->m, it2_min->m, row, fullSide);
                            if(it1_min->is_fixed && it2_min->is_fixed){
                                tmp = matSum(it1_min->m, it2_min->m);
                            }else{
                                tmp = matSum1(row, it1_min->m, it2_min->m, col);
                            }
                            if (it1_min->is_tmp) matDestroy(it1_min->m);
                            if (it2_min->is_tmp) matDestroy(it2_min->m);

                            rl.insert(it1_min, data_type{tmp, true, true});
                            rl.erase(it1_min);
                            rl.erase(it2_min);
                        }
                        it1_min = it2_min = rl.begin();
                        ++it2_min;
                        //tmp = matSum1(it1_min->m, it2_min->m, row, fullSide);
                        if(it1_min->is_fixed && it2_min->is_fixed){
                            tmp = matSum(it1_min->m, it2_min->m);
                        }else{
                            tmp = matSum1(row, it1_min->m, it2_min->m, col);
                        }
                        if(it1_min->is_tmp) matDestroy(it1_min->m);
                        if(it2_min->is_tmp) matDestroy(it2_min->m);
                        res.insert(res.begin(), data_type{tmp, true, true});
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, STAR, ll);
                    if(!skip_closure){
                        matrix tmp = matClos1(row, ll.front().m, 0, col);
                        // std::cout << "STAR : " << ll.front() << std::endl;
                        if(ll.front().is_tmp) matDestroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, true, true});
                    }else{
                        res = std::move(ll);
                    }

                    break;
                }
                case PLUS:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, PLUS, ll);
                    if(!skip_closure){
                        matrix tmp = matClos1(row, ll.front().m, 1, col);
                        if(ll.front().is_tmp) matDestroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, true, true});
                    }else{
                        res = std::move(ll);
                    }
                    break;
                }
                case QUESTION:
                {
                    list_type ll;
                    traversal_row_col_fixed(rpqTree, node->e1, QUESTION, row, col, ll);
                    matrix Id = matId(std::max(ll.front().m->height,ll.front().m->width));
                    matrix tmp = matSum1(row, ll.front().m, Id, col);
                    matDestroy(Id);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true, true});
                    break;
                }
            }
        }

        std::string file_name(const uint i, const std::string &index){
            std::stringstream ss;
            ss << std::setw(4) << std::setfill('0') << i;
            return index + "/" + ss.str() + ".mat";
        }

    public:

        //TODO: should be const
        std::unordered_map<std::string, uint64_t> &map_SO = m_map_SO;
        std::unordered_map<std::string, uint64_t> &map_P = m_map_P;

        explicit solver(const std::string &dataset, const std::string &index){

            m_matrices.resize(SIZE);
            uint64_t space = 0;
            uint i;
            FILE *f;
            printf ("Reading 5420 matrices..."); fflush(stdout);
            for (i=1;i<SIZE;i++) {
                std::string file = file_name(i, index);
                f = fopen(file.c_str(), "r");
                m_matrices[i] = matLoad(f);
                fclose(f);
                space += matSpace(m_matrices[i]);
            }
            printf (" done... %li total words (%0.2f bpt)\n",space,space*(w/8)/(float)N);


            std::ifstream ifs_SO(dataset + ".SO", std::ifstream::in);
            std::ifstream ifs_P(dataset + ".P", std::ifstream::in);
            //std::ifstream ifs_q(argv[2], std::ifstream::in);
            //std::ofstream ofs_SO("SO", std::ofstream::out);


            uint64_t id;
            std::string s_aux, data;

            std::cout << "Reading mapping..." << std::flush;
            while (std::getline(ifs_SO, data)) {
                auto space = data.find(' ');
                id = std::stoull(data.substr(0, space));
                s_aux = data.substr(space + 1);
                m_map_SO[s_aux] = id;
            }

            while (std::getline(ifs_P, data)) {
                auto space = data.find(' ');
                id = std::stoull(data.substr(0, space));
                s_aux = data.substr(space + 1);
                m_map_P[s_aux] = id;
            }
            std::cout << " done." << std::endl;

        }

        matrix solve_var_to_var(std::string &query){
            list_type res;
            RpqTree rpqTree(query, map_P, SIZE);
            traversal(&rpqTree, rpqTree.root(), ROOT, res);
            return res.front().m;
        }

        matrix solve_con_to_var(std::string &query, int s_id){
            list_type res;
            RpqTree rpqTree(query, map_P, SIZE);
            traversal_row_fixed(&rpqTree, rpqTree.root(), ROOT, s_id, res);
            return res.front().m;
        }

        matrix solve_var_to_con(std::string &query, int o_id){
            list_type res;
            RpqTree rpqTree(query, map_P, SIZE);
            traversal_col_fixed(&rpqTree, rpqTree.root(), ROOT, o_id, res);
            return res.front().m;
        }

        matrix solve_con_to_con(std::string &query, int s_id, int o_id){
            list_type res;
            RpqTree rpqTree(query, map_P, SIZE);
            traversal_row_col_fixed(&rpqTree, rpqTree.root(), ROOT, s_id, o_id, res);
            return res.front().m;
        }

    };


}

#endif //MATRIX_RPQ_RPQ_SOLVER_HPP