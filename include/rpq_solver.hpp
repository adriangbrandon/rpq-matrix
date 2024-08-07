//
// Created by Adrián on 3/5/23.
//

#ifndef MATRIX_RPQ_RPQ_SOLVER_HPP
#define MATRIX_RPQ_RPQ_SOLVER_HPP

#define N 958844164
#define SIZE 5420 // 1 to 5419
#define V 296008192 // 1 to...

#ifndef w
#define w (8*sizeof(uint64_t))
#endif

#define VERBOSE 0


#include "RpqTree.hpp"
#include <list>
#include <vector>

namespace rpq {

    template <class wrapper>
    class solver {

    public:

        typedef typename wrapper::matrix_type matrix;
        typedef typename wrapper::s_matrix_type s_matrix;
        typedef struct {
            matrix m;
            bool is_transposed;
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


        inline matrix get_matrix(s_matrix &a, matrix m, bool is_tranposed){
            if(!is_tranposed) return m;
            a = wrapper::transpose(m);
            return &a;
        }

        void traversal(RpqTree* rpqTree, Tree* node, int parentType, list_type &res,
                       bool skip_closure = false){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    matrix a = (pred > SIZE) ? m_matrices[pred-SIZE] :  m_matrices[pred];
                    res.insert(res.begin(), data_type{a, (pred > SIZE),
                                                      false, false});
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
#if VERBOSE
                        std::cout << "[CONC]: skip" << std::endl;
#endif
                    }else{
                        it_type it1, it2, it1_min, it2_min;
                        matrix tmp;
                        while(rl.size() > 2){//Check size if there are more than three elements
                            it1 = it2 = rl.begin();
                            ++it2;
                            uint64_t min = UINT64_MAX, weight;
                            while(it2 != rl.end()){
                                weight = wrapper::space(it1->m) + wrapper::space(it2->m);
                                if(min > weight){
                                    it1_min = it1;
                                    it2_min = it2;
                                    min = weight;
                                }
                                ++it1; ++it2;
                            }
                            matrix A, B;
                            s_matrix sA, sB;
                            A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                            B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                            tmp = wrapper::mult(A, B);
                            if(it1_min->is_tmp) wrapper::destroy(it1_min->m);
                            if(it2_min->is_tmp) wrapper::destroy(it2_min->m);

                            rl.insert(it1_min, data_type{tmp, false, true, false});
                            rl.erase(it1_min);
                            rl.erase(it2_min);
                        }
                        it1_min = it2_min = rl.begin();
                        ++it2_min;
                        matrix A, B;
                        s_matrix sA, sB;
                        A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                        B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                        tmp = wrapper::mult(A, B);
                        if(it1_min->is_tmp) wrapper::destroy(it1_min->m);
                        if(it2_min->is_tmp) wrapper::destroy(it2_min->m);
                        res.insert(res.begin(), data_type{tmp, false, true, false});
#if VERBOSE
                        std::cout << "[CONC]: wrapper::mult" << std::endl;
#endif
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
#if VERBOSE
                        std::cout << "[OR]: skip" << std::endl;
#endif
                    }else{
                        it_type it1, it2, it1_min, it2_min;
                        matrix tmp;
                        while(rl.size() > 2){//Check size if there are more than three elements
                            it1 = it2 = rl.begin();
                            uint64_t min = UINT64_MAX, weight;
                            while(it1 != rl.end()){
                                weight = wrapper::space(it1->m);
                                if(min > weight){
                                    it1_min = it1;
                                    min = weight;
                                }
                                ++it1;
                            }
                            min = UINT64_MAX;
                            while(it2 != rl.end()){
                                if(it2 != it1_min ){
                                    weight = wrapper::space(it2->m);
                                    if(min > weight){
                                        it2_min = it2;
                                        min = weight;
                                    }
                                }
                                ++it2;
                            }
                            matrix A, B;
                            s_matrix sA, sB;
                            A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                            B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                            tmp = wrapper::sum(A, B);
                            if(it1_min->is_tmp) wrapper::destroy(it1_min->m);
                            if(it2_min->is_tmp) wrapper::destroy(it2_min->m);

                            rl.insert(it1_min, data_type{tmp, false, true, false});
                            rl.erase(it1_min);
                            rl.erase(it2_min);
                        }

                        it1_min = it2_min = rl.begin();
                        ++it2_min;
                        matrix A, B;
                        s_matrix sA, sB;
                        A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                        B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                        tmp = wrapper::sum(A, B);
                        if(it1_min->is_tmp) wrapper::destroy(it1_min->m);
                        if(it2_min->is_tmp) wrapper::destroy(it2_min->m);
                        res.insert(res.begin(), data_type{tmp, false,true, false});
#if VERBOSE
                        std::cout << "[OR]: wrapper::sum" << std::endl;
#endif
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, STAR, ll);
                    if(!skip_closure){
                        matrix A;
                        s_matrix sA;
                        A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                        matrix tmp = wrapper::clos(A, 0);
                        if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                        // std::cout << "STAR : " << ll.front() << std::endl;
                        res.insert(res.begin(), data_type{tmp, false, true, false});
#if VERBOSE
                        std::cout << "[STAR]: matClos" << std::endl;
#endif
                    }else{
                        res = std::move(ll);
#if VERBOSE
                        std::cout << "[STAR]: skip" << std::endl;
#endif
                    }
                    break;
                }
                case PLUS:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, PLUS, ll);
                    if(!skip_closure){
                        matrix A;
                        s_matrix sA;
                        A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                        matrix tmp = wrapper::clos(A, 1);
                        if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, false,true, false});
                    }else{
                        res = std::move(ll);
                    }
                    break;
                }
                case QUESTION:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, QUESTION, ll);
                    matrix A;
                    s_matrix sA;
                    A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                    matrix Id = wrapper::id(std::max(A->height, A->width));
                    matrix tmp = wrapper::sum(A, Id);
                    wrapper::destroy(Id);
                    if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, false,true, false});
                    break;
                }
            }
        }

        void traversal_col_fixed(RpqTree* rpqTree, Tree* node, int parentType, int col, list_type &res,
                                 bool skip_closure = false){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    if(parentType != ROOT){
                        matrix a = (pred > SIZE) ? m_matrices[pred-SIZE] :  m_matrices[pred];
                        res.insert(res.begin(), data_type{a, (pred > SIZE), false, false});
                    }else {
                        matrix a = (pred > SIZE) ? m_matrices[pred-SIZE] :  m_matrices[pred];
                        matrix A;
                        s_matrix sA;
                        A = get_matrix(sA, a, (pred > SIZE));
                        matrix e = wrapper::empty(A->height, A->width);
                        matrix m = wrapper::sum1(wrapper::full_side, A, e, col);
                        wrapper::destroy(e);
                        res.insert(res.begin(), data_type{m, false,true, true});
                    }
                    break;
                }
                case CONC:
                {
                    list_type ll, rl;
                    if(node->e1->type == STAR || node->e1->type == PLUS){ //Special case!
                        traversal_col_fixed(rpqTree, node->e2, CONC, col, rl);
                        traversal(rpqTree, node->e1, CONC, ll, true);
                        matrix tmp;
                        uint pos = (node->e1->type == STAR) ? 0 : 1;
                        if(rl.size() > 1){ //product of the right part
                            auto it2 = rl.rbegin(); //last element
                            auto it1 = rl.rbegin(); //last element
                            ++it2; //last element -1
                            matrix aux;
                            if(it1->is_fixed){
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it2->m, it2->is_transposed);
                                B = get_matrix(sB, it1->m, it1->is_transposed);
                                tmp = wrapper::mult(A, B);
                            }else{
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it2->m, it2->is_transposed);
                                B = get_matrix(sB, it1->m, it1->is_transposed);
                                tmp = wrapper::mult1(wrapper::full_side, A,B, col);
                            }
                            if(it1->is_tmp) wrapper::destroy(it1->m);
                            if(it2->is_tmp) wrapper::destroy(it2->m);
                            while(++it2 != rl.rend()){
                                matrix A;
                                s_matrix sA;
                                A = get_matrix(sA, it2->m, it2->is_transposed);
                                aux = wrapper::mult(A, tmp); //tmp is already fixed
                                if(it2->is_tmp) wrapper::destroy(it2->m);
                                wrapper::destroy(tmp);
                                tmp = aux;
                            }
                            matrix A;
                            s_matrix sA;
                            A = get_matrix(sA, ll.front().m,  ll.front().is_transposed);
                            aux = wrapper::clos_mult1(wrapper::full_side, A, pos, tmp, col);
                            wrapper::destroy(tmp);
                            if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                            res.insert(res.begin(), data_type{aux, false, true, true});
                        }else{
                            matrix A, B;
                            s_matrix sA, sB;
                            A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                            B = get_matrix(sB, rl.front().m, rl.front().is_transposed);
                            tmp = wrapper::clos_mult1(wrapper::full_side, A, pos, B, col);
                            if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                            if(rl.front().is_tmp) wrapper::destroy(rl.front().m);
                            res.insert(res.begin(), data_type{tmp, false, true, true});
                        }
                    }else{
                        traversal_col_fixed(rpqTree, node->e2, CONC, col, rl);
                        traversal(rpqTree, node->e1, CONC, ll);
                        rl.splice(rl.begin(), ll);
                        if (parentType == CONC){
                            res = std::move(rl);
                        }else{
                            auto it2 = rl.rbegin(); //last element
                            auto it1 = rl.rbegin(); //last element
                            ++it2; //last element -1
                            matrix tmp, aux;
                            if(it1->is_fixed){
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it2->m, it2->is_transposed);
                                B = get_matrix(sB, it1->m, it1->is_transposed);
                                tmp = wrapper::mult(A, B);
                            }else{
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it2->m, it2->is_transposed);
                                B = get_matrix(sB, it1->m, it1->is_transposed);
                                tmp = wrapper::mult1(wrapper::full_side, A, B, col);
                            }
                            if(it1->is_tmp) wrapper::destroy(it1->m);
                            if(it2->is_tmp) wrapper::destroy(it2->m);
                            //TODO: original
                            /*while(++it2 != rl.rend()){
                                matrix A;
                                s_matrix sA;
                                A = get_matrix(sA, it2->m, it2->is_transposed);
                                aux = wrapper::mult(A, tmp); //tmp is already fixed
                                if(it2->is_tmp) wrapper::destroy(it2->m);
                                wrapper::destroy(tmp);
                                tmp = aux;
                            }
                            res.insert(res.begin(), data_type{tmp, false, true, true});
                            */

                            if(rl.size() == 2) {
                                res.insert(res.begin(), data_type{tmp, false, true, true});
                            }else {
                                rl.insert(it2, data_type{tmp, false, true, true});
                                rl.erase(it2);
                                rl.erase(it1);
                                it_type it1_min, it2_min;
                                while(rl.size() > 2){//Check size if there are more than three elements
                                    it1 = it2 = rl.begin();
                                    ++it2;
                                    uint64_t min = UINT64_MAX, weight;
                                    while(it2 != rl.end()){
                                        weight = wrapper::space(it1->m) + wrapper::space(it2->m);
                                        if(min > weight){
                                            it1_min = it1;
                                            it2_min = it2;
                                            min = weight;
                                        }
                                        ++it1; ++it2;
                                    }
                                    matrix A, B;
                                    s_matrix sA, sB;
                                    A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                                    B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                                    tmp = wrapper::mult(A, B);
                                    if(it1_min->is_tmp) wrapper::destroy(it1_min->m);
                                    if(it2_min->is_tmp) wrapper::destroy(it2_min->m);

                                    rl.insert(it1_min, data_type{tmp, false, true, false});
                                    rl.erase(it1_min);
                                    rl.erase(it2_min);
                                }
                                it1_min = it2_min = rl.begin();
                                ++it2_min;
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                                B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                                tmp = wrapper::mult(A, B);
                                if(it1_min->is_tmp) wrapper::destroy(it1_min->m);
                                if(it2_min->is_tmp) wrapper::destroy(it2_min->m);

                                res.insert(res.begin(), data_type{tmp, false, true, true});
                            }
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
#if VERBOSE
                        std::cout << "[OR]: col fixed -> skip" << std::endl;
#endif
                    }else{
                        matrix tmp;
                        it_type it1, it2, it1_min, it2_min;
                        while(rl.size() > 2) {//Check if there a
                            it1 = it2 = rl.begin();
                            uint64_t min = UINT64_MAX, weight;
                            while (it1 != rl.end()) {
                                weight = wrapper::space(it1->m);
                                if (min > weight) {
                                    it1_min = it1;
                                    min = weight;
                                }
                                ++it1; //++it2;
                            }
                            min = UINT64_MAX;
                            while (it2 != rl.end()) {
                                if (it2 != it1_min) {
                                    weight = wrapper::space(it2->m);
                                    if (min > weight) {
                                        it2_min = it2;
                                        min = weight;
                                    }
                                }
                                ++it2;
                            }
                            //tmp = matSum1(it1_min->m, it2_min->m, wrapper::full_side, col);
                            if(it1_min->is_fixed && it2_min->is_fixed){
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                                B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                                tmp = wrapper::sum(A, B);
                            }else{
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                                B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                                tmp = wrapper::sum1(wrapper::full_side, A, B, col);
                            }
                            if (it1_min->is_tmp) wrapper::destroy(it1_min->m);
                            if (it2_min->is_tmp) wrapper::destroy(it2_min->m);

                            rl.insert(it1_min, data_type{tmp, false,true, true});
                            rl.erase(it1_min);
                            rl.erase(it2_min);
                        }
                        it1_min = it2_min = rl.begin();
                        ++it2_min;
                        if(it1_min->is_fixed && it2_min->is_fixed){
                            matrix A, B;
                            s_matrix sA, sB;
                            A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                            B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                            tmp = wrapper::sum(A, B);
                        }else{
                            matrix A, B;
                            s_matrix sA, sB;
                            A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                            B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                            tmp = wrapper::sum1(wrapper::full_side, A, B, col);
                        }
                        if(it1_min->is_tmp) wrapper::destroy(it1_min->m);
                        if(it2_min->is_tmp) wrapper::destroy(it2_min->m);
                        res.insert(res.begin(), data_type{tmp, false,true, true});
#if VERBOSE
                        std::cout << "[OR]: col fixed -> matSums" << std::endl;
#endif
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, STAR, ll);
                    if(!skip_closure){
#if VERBOSE
                        std::cout << "[STAR]: col fixed -> matClos1" << std::endl;
#endif
                        matrix A;
                        s_matrix sA;
                        A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                        matrix tmp = wrapper::clos1(wrapper::full_side, A, 0, col);
                        // std::cout << "STAR : " << ll.front() << std::endl;
                        if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, false, true, true});
                    }else{
                        res = std::move(ll);
#if VERBOSE
                        std::cout << "[STAR]: col fixed -> skip" << std::endl;
#endif
                    }

                    break;
                }
                case PLUS:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, PLUS, ll);
                    if(!skip_closure){
#if VERBOSE
                        std::cout << "[PLUS]: col fixed -> matClos1" << std::endl;
#endif
                        matrix A;
                        s_matrix sA;
                        A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                        matrix tmp = wrapper::clos1(wrapper::full_side, A, 1, col);
                        if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, false,true, true});
                    }else{
#if VERBOSE
                        std::cout << "[PLUS]: col fixed -> skip" << std::endl;
#endif
                        res = std::move(ll);
                    }
                    break;
                }
                case QUESTION:
                {

                    list_type ll;
                    traversal_col_fixed(rpqTree, node->e1, QUESTION, col, ll);
                    matrix A;
                    s_matrix sA;
                    A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                    matrix Id = wrapper::id1(std::max(A->height,A->width), col);
                    matrix tmp = wrapper::sum1(wrapper::full_side, A, Id, col);
                    wrapper::destroy(Id);
                    if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, false, true, true});
                    break;
                }
            }
        }

        void traversal_row_fixed(RpqTree* rpqTree, Tree* node, int parentType, int row, list_type &res,
                                 bool skip_closure = false){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    if(parentType != ROOT){
                        matrix a = (pred > SIZE) ? m_matrices[pred-SIZE] :  m_matrices[pred];
                        res.insert(res.begin(), data_type{a, (pred > SIZE), false, false});
                    }else {
                        matrix A;
                        s_matrix sA;
                        matrix a = (pred > SIZE) ? m_matrices[pred-SIZE] :  m_matrices[pred];
                        A = get_matrix(sA, a, (pred > SIZE));
                        matrix e = wrapper::empty(A->height, A->width);
                        matrix m = wrapper::sum1(row, A, e, wrapper::full_side);
                        wrapper::destroy(e);
                        res.insert(res.begin(), data_type{m, false, true, true});
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
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1->m, it1->is_transposed);
                                B = get_matrix(sB, it2->m, it2->is_transposed);
                                tmp = wrapper::mult(A, B);
                            }else{
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1->m, it1->is_transposed);
                                B = get_matrix(sB, it2->m, it2->is_transposed);
                                tmp = wrapper::mult1(row, A, B, wrapper::full_side);
                            }
                            if(it1->is_tmp) wrapper::destroy(it1->m);
                            if(it2->is_tmp) wrapper::destroy(it2->m);
                            while(++it2 != ll.end()){
                                matrix A;
                                s_matrix sA;
                                A = get_matrix(sA, it2->m, it2->is_transposed);
                                aux = wrapper::mult(A, tmp); //tmp is already fixed
                                if(it2->is_tmp) wrapper::destroy(it2->m);
                                wrapper::destroy(tmp);
                                tmp = aux;
                            }
                            aux = wrapper::mult_clos1(row, tmp, rl.front().m, pos, wrapper::full_side);
                            wrapper::destroy(tmp);
                            if(rl.front().is_tmp) wrapper::destroy(rl.front().m);
                            res.insert(res.begin(), data_type{aux, false, true, true});
                        }else{
                            matrix A, B;
                            s_matrix sA, sB;
                            A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                            B = get_matrix(sB, rl.front().m, rl.front().is_transposed);
                            tmp = wrapper::mult_clos1(row, A, B, pos, wrapper::full_side);
                            if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                            if(rl.front().is_tmp) wrapper::destroy(rl.front().m);
                            res.insert(res.begin(), data_type{tmp, false,true, true});
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
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1->m, it1->is_transposed);
                                B = get_matrix(sB, it2->m, it2->is_transposed);
                                tmp = wrapper::mult(A, B);
                            }else{
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1->m, it1->is_transposed);
                                B = get_matrix(sB, it2->m, it2->is_transposed);
                                tmp = wrapper::mult1(row, A,B, wrapper::full_side);
                            }
                            if(it1->is_tmp) wrapper::destroy(it1->m);
                            if(it2->is_tmp) wrapper::destroy(it2->m);
                            //TODO: original
                            /*while(++it2 != rl.end()){
                                matrix A;
                                s_matrix sA;
                                A = get_matrix(sA, it2->m, it2->is_transposed);
                                aux = wrapper::mult(tmp, A);
                                if(it2->is_tmp) wrapper::destroy(it2->m);
                                wrapper::destroy(tmp);
                                tmp = aux;
                            }
                            res.insert(res.begin(), data_type{tmp, false,true, true});
                            */
                            if(rl.size() == 2) {
                                res.insert(res.begin(), data_type{tmp, false, true, true});
                            }else {
                                rl.insert(it1, data_type{tmp, false, true, true});
                                rl.erase(it1);
                                rl.erase(it2);
                                it_type it1_min, it2_min;
                                while(rl.size() > 2){//Check size if there are more than three elements
                                    it1 = it2 = rl.begin();
                                    ++it2;
                                    uint64_t min = UINT64_MAX, weight;
                                    while(it2 != rl.end()){
                                        weight = wrapper::space(it1->m) + wrapper::space(it2->m);
                                        if(min > weight){
                                            it1_min = it1;
                                            it2_min = it2;
                                            min = weight;
                                        }
                                        ++it1; ++it2;
                                    }
                                    matrix A, B;
                                    s_matrix sA, sB;
                                    A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                                    B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                                    tmp = wrapper::mult(A, B);
                                    if(it1_min->is_tmp) wrapper::destroy(it1_min->m);
                                    if(it2_min->is_tmp) wrapper::destroy(it2_min->m);

                                    rl.insert(it1_min, data_type{tmp, false, true, false});
                                    rl.erase(it1_min);
                                    rl.erase(it2_min);
                                }
                                it1_min = it2_min = rl.begin();
                                ++it2_min;
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                                B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                                tmp = wrapper::mult(A, B);
                                if(it1_min->is_tmp) wrapper::destroy(it1_min->m);
                                if(it2_min->is_tmp) wrapper::destroy(it2_min->m);

                                res.insert(res.begin(), data_type{tmp, false, true, true});
                            }
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
                                weight = wrapper::space(it1->m);
                                if (min > weight) {
                                    it1_min = it1;
                                    min = weight;
                                }
                                ++it1; //++it2;
                            }
                            min = UINT64_MAX;
                            while (it2 != rl.end()) {
                                if (it2 != it1_min) {
                                    weight = wrapper::space(it2->m);
                                    if (min > weight) {
                                        it2_min = it2;
                                        min = weight;
                                    }
                                }
                                ++it2;
                            }
                            //tmp = wrapper::sum1(it1_min->m, it2_min->m, row, wrapper::full_side);
                            if(it1_min->is_fixed && it2_min->is_fixed){
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                                B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                                tmp = wrapper::sum(A, B);
                            }else{
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                                B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                                tmp = wrapper::sum1(row, A, B, wrapper::full_side);
                            }
                            if (it1_min->is_tmp) wrapper::destroy(it1_min->m);
                            if (it2_min->is_tmp) wrapper::destroy(it2_min->m);

                            rl.insert(it1_min, data_type{tmp, false, true, true});
                            rl.erase(it1_min);
                            rl.erase(it2_min);
                        }
                        it1_min = it2_min = rl.begin();
                        ++it2_min;
                        //tmp = matSum1(it1_min->m, it2_min->m, row, wrapper::full_side);
                        if(it1_min->is_fixed && it2_min->is_fixed){
                            matrix A, B;
                            s_matrix sA, sB;
                            A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                            B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                            tmp = wrapper::sum(A, B);
                        }else{
                            matrix A, B;
                            s_matrix sA, sB;
                            A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                            B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                            tmp = wrapper::sum1(row, A, B, wrapper::full_side);
                        }
                        if(it1_min->is_tmp) wrapper::destroy(it1_min->m);
                        if(it2_min->is_tmp) wrapper::destroy(it2_min->m);
                        res.insert(res.begin(), data_type{tmp, false, true, true});
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, STAR, ll);
                    if(!skip_closure){
                        matrix A;
                        s_matrix sA;
                        A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                        matrix tmp = wrapper::clos1(row, A, 0, wrapper::full_side);
                        // std::cout << "STAR : " << ll.front() << std::endl;
                        if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, false, true, true});
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
                        matrix A;
                        s_matrix sA;
                        A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                        matrix tmp = wrapper::clos1(row, A, 1, wrapper::full_side);
                        if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, false, true, true});
                    }else{
                        res = std::move(ll);
                    }
                    break;
                }
                case QUESTION:
                {
                    list_type ll;
                    traversal_row_fixed(rpqTree, node->e1, QUESTION, row, ll);
                    matrix A;
                    s_matrix sA;
                    A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                    matrix Id = wrapper::id1(std::max(A->height,A->width), row);
                    matrix tmp = wrapper::sum1(row, ll.front().m, Id, wrapper::full_side);
                    wrapper::destroy(Id);
                    if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, false, true, true});
                    break;
                }
            }
        }

        void traversal_row_col_fixed(RpqTree* rpqTree, Tree* node, int parentType, int row, int col, list_type &res,
                                     bool skip_closure = false){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    if(parentType != ROOT){
                        matrix a = (pred > SIZE) ? m_matrices[pred-SIZE] :  m_matrices[pred];
                        res.insert(res.begin(), data_type{a, (pred > SIZE), false, false});
                    }else {
                        matrix A;
                        s_matrix sA;
                        matrix a = (pred > SIZE) ? m_matrices[pred-SIZE] :  m_matrices[pred];
                        A = get_matrix(sA, a, (pred > SIZE));
                        matrix e = wrapper::empty(A->height, A->width);
                        matrix m = wrapper::sum1(row, A, e, col);
                        wrapper::destroy(e);
                        res.insert(res.begin(), data_type{m, false, true, true});
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
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1->m, it1->is_transposed);
                                B = get_matrix(sB, it2->m, it2->is_transposed);
                                tmp = wrapper::mult(A, B);
                            }else{
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1->m, it1->is_transposed);
                                B = get_matrix(sB, it2->m, it2->is_transposed);
                                tmp = wrapper::mult1(row, A, B, wrapper::full_side);
                            }
                            if(it1->is_tmp) wrapper::destroy(it1->m);
                            if(it2->is_tmp) wrapper::destroy(it2->m);
                            while(++it2 != ll.end()){
                                matrix A;
                                s_matrix sA;
                                A = get_matrix(sA, it2->m, it2->is_transposed);
                                aux = wrapper::mult(A, tmp); //tmp is already fixed
                                if(it2->is_tmp) wrapper::destroy(it2->m);
                                wrapper::destroy(tmp);
                                tmp = aux;
                            }
                            aux = wrapper::mult_clos1(row, tmp, rl.front().m, pos, col);
                            wrapper::destroy(tmp);
                            if(rl.front().is_tmp) wrapper::destroy(rl.front().m);
                            res.insert(res.begin(), data_type{aux, false, true, true});
                        }else{
                            matrix A, B;
                            s_matrix sA, sB;
                            A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                            B = get_matrix(sB, rl.front().m, rl.front().is_transposed);
                            tmp = wrapper::mult_clos1(row, A, B, pos, col);
                            if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                            if(rl.front().is_tmp) wrapper::destroy(rl.front().m);
                            res.insert(res.begin(), data_type{tmp, false,true, true});
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
                            uint64_t cnt_mult = 2;
                            if(it1->is_fixed){
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1->m, it1->is_transposed);
                                B = get_matrix(sB, it2->m, it2->is_transposed);
                                if(rl.size() == cnt_mult){
                                    tmp = wrapper::mult1(row, A, B, col); //Final
                                }else {
                                    tmp = wrapper::mult(A, B);
                                }
                            }else{
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1->m, it1->is_transposed);
                                B = get_matrix(sB, it2->m, it2->is_transposed);
                                if(rl.size() == cnt_mult){
                                    tmp = wrapper::mult1(row, A, B, col); //Final
                                }else{
                                    tmp = wrapper::mult1(row, A,B, wrapper::full_side);
                                }
                            }
                            if(it1->is_tmp) wrapper::destroy(it1->m);
                            if(it2->is_tmp) wrapper::destroy(it2->m);
                            while(++it2 != rl.end()){
                                ++cnt_mult;
                                matrix A;
                                s_matrix sA;
                                A = get_matrix(sA, it2->m, it2->is_transposed);
                                if(rl.size() == cnt_mult){
                                    aux = wrapper::mult1(row, tmp, A, col); //Final
                                }else{
                                    aux = wrapper::mult(tmp, A);
                                }
                                if(it2->is_tmp) wrapper::destroy(it2->m);
                                wrapper::destroy(tmp);
                                tmp = aux;
                            }
                            res.insert(res.begin(), data_type{tmp, false,true, true});
                        }
                    }
                    break;
                }
                case OOR:
                {
                    list_type ll, rl;
                    traversal_row_col_fixed(rpqTree, node->e1, OOR, row, col,ll);
                    traversal_row_col_fixed(rpqTree, node->e2, OOR, row, col,rl);
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
                                weight = wrapper::space(it1->m);
                                if (min > weight) {
                                    it1_min = it1;
                                    min = weight;
                                }
                                ++it1; //++it2;
                            }
                            min = UINT64_MAX;
                            while (it2 != rl.end()) {
                                if (it2 != it1_min) {
                                    weight = wrapper::space(it2->m);
                                    if (min > weight) {
                                        it2_min = it2;
                                        min = weight;
                                    }
                                }
                                ++it2;
                            }
                            //tmp = wrapper::sum1(it1_min->m, it2_min->m, row, wrapper::full_side);
                            if(it1_min->is_fixed && it2_min->is_fixed){
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                                B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                                tmp = wrapper::sum(A, B);
                            }else{
                                matrix A, B;
                                s_matrix sA, sB;
                                A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                                B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                                tmp = wrapper::sum1(row, A, B, wrapper::full_side);
                            }
                            if (it1_min->is_tmp) wrapper::destroy(it1_min->m);
                            if (it2_min->is_tmp) wrapper::destroy(it2_min->m);

                            rl.insert(it1_min, data_type{tmp, false, true, true});
                            rl.erase(it1_min);
                            rl.erase(it2_min);
                        }
                        it1_min = it2_min = rl.begin();
                        ++it2_min;
                        //tmp = wrapper::sum1(it1_min->m, it2_min->m, row, wrapper::full_side);
                        if(it1_min->is_fixed && it2_min->is_fixed){
                            matrix A, B;
                            s_matrix sA, sB;
                            A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                            B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                            tmp = wrapper::sum(A, B);
                        }else{
                            matrix A, B;
                            s_matrix sA, sB;
                            A = get_matrix(sA, it1_min->m, it1_min->is_transposed);
                            B = get_matrix(sB, it2_min->m, it2_min->is_transposed);
                            tmp = wrapper::sum1(row, A, B, col);
                        }
                        if(it1_min->is_tmp) wrapper::destroy(it1_min->m);
                        if(it2_min->is_tmp) wrapper::destroy(it2_min->m);
                        res.insert(res.begin(), data_type{tmp, false, true, true});
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, STAR, ll);
                    if(!skip_closure){
                        matrix A;
                        s_matrix sA;
                        A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                        matrix tmp = wrapper::clos1(row, A, 0, col);
                        // std::cout << "STAR : " << ll.front() << std::endl;
                        if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, false, true, true});
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
                        matrix A;
                        s_matrix sA;
                        A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                        matrix tmp = wrapper::clos1(row, A, 1, col);
                        if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                        res.insert(res.begin(), data_type{tmp, false, true, true});
                    }else{
                        res = std::move(ll);
                    }
                    break;
                }
                case QUESTION:
                {
                    list_type ll;
                    traversal_row_col_fixed(rpqTree, node->e1, QUESTION, row, col,ll);
                    matrix A;
                    s_matrix sA;
                    A = get_matrix(sA, ll.front().m, ll.front().is_transposed);
                    matrix Id = wrapper::id1(std::max(A->height,A->width), col);
                    matrix tmp = wrapper::sum1(row, ll.front().m, Id, col);
                    wrapper::destroy(Id);
                    if(ll.front().is_tmp) wrapper::destroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, false, true, true});
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

        explicit solver(const std::string &dataset, const std::string &index,
                        const uint n_preds, const uint n_triples){

            m_matrices.resize(n_preds + 1);
            uint64_t space = 0;
            uint i;
            FILE *f;
            printf ("Reading %i matrices...", n_preds); fflush(stdout);
            for (i=1;i<m_matrices.size();i++) {
                std::string file = file_name(i, index);
                f = fopen(file.c_str(), "r");
                m_matrices[i] = wrapper::load(f);
                fclose(f);
                space += wrapper::space(m_matrices[i]);
            }
            printf (" done... %li total words (%0.2f bpt)\n",space,space*(w/8)/(float)n_triples);

            std::ifstream ifs_SO(dataset + ".SO", std::ifstream::in);
            std::ifstream ifs_P(dataset + ".P", std::ifstream::in);
            //std::ifstream ifs_q(argv[2], std::ifstream::in);
            //std::ofstream ofs_SO("SO", std::ofstream::out);


            uint64_t id;
            std::string s_aux, data;

            std::cout << "Reading mapping..." << std::flush;
            while (std::getline(ifs_SO, data)) {
                space = data.find(' ');
                id = std::stoull(data.substr(0, space));
                s_aux = data.substr(space + 1);
                m_map_SO[s_aux] = id;
            }

            while (std::getline(ifs_P, data)) {
                space = data.find(' ');
                id = std::stoull(data.substr(0, space));
                s_aux = data.substr(space + 1);
                m_map_P[s_aux] = id;
            }
            std::cout << " done." << std::endl;

            /*double_t bits = std::ceil(std::log2(N)) * (SIZE-1);
            for(uint64_t j = 1; j < m_matrices.size(); ++j){
                bits = bits + std::ceil(std::log2(m_matrices[j]->elems)) * m_matrices[j]->elems;
            }
            double_t bytes = bits / 8;
            double_t bpt = bytes / N;
            std::cout << "====== VP =====" << std::endl;
            std::cout << bits << " (bits)" << std::endl;
            std::cout << bytes << " (bytes)" << std::endl;
            std::cout << bpt << " (bpt)" << std::endl;*/

        }

        data_type solve_var_to_var(std::string &query, bool &rem){
            list_type res;
            RpqTree rpqTree(query, map_P, m_matrices.size());
            traversal(&rpqTree, rpqTree.root(), ROOT, res);
            return res.front();
        }

        data_type solve_con_to_var(std::string &query, int s_id, bool &rem){
            list_type res;
            RpqTree rpqTree(query, map_P, m_matrices.size());
            traversal_row_fixed(&rpqTree, rpqTree.root(), ROOT, s_id, res);
            return res.front();
        }

        data_type solve_var_to_con(std::string &query, int o_id, bool &rem){
            list_type res;
            RpqTree rpqTree(query, map_P, m_matrices.size());
            traversal_col_fixed(&rpqTree, rpqTree.root(), ROOT, o_id, res);
            return res.front();
        }

        data_type solve_con_to_con(std::string &query, int s_id, int o_id, bool &rem){
            list_type res;
            RpqTree rpqTree(query, map_P, m_matrices.size());
            traversal_row_col_fixed(&rpqTree, rpqTree.root(), ROOT, s_id, o_id, res);
            return res.front();
        }

    };


}

#endif //MATRIX_RPQ_RPQ_SOLVER_HPP