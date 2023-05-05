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
        } data_type;
        typedef std::list<data_type> list_type;
        typedef typename list_type::iterator it_type;
        typedef typename list_type::reverse_iterator rev_it_type;

    private:

        std::vector<matrix> m_matrices;
        //std::vector<matrix> m_tmp_matrices;
        std::unordered_map<std::string, uint64_t> m_map_SO;
        std::unordered_map<std::string, uint64_t> m_map_P;

        void traversal(RpqTree* rpqTree, Tree* node, int parentType, list_type &res){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    if(pred > SIZE){
                        std::cout << "Neg." << std::endl;
                    }else {
                        res.insert(res.begin(), data_type{m_matrices[pred], false});
                    }
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
                        while(rl.size() > 1){//Check if there is only one element
                            it1 = it2 = rl.begin();
                            ++it2;
                            uint64_t min = UINT64_MAX, weight;
                            while(it2 != rl.end()){
                                weight = std::min(matSpace(it1->m), matSpace(it2->m));
                                if(min > weight){
                                    it1_min = it1;
                                    it2_min = it2;
                                    min = weight;
                                }
                                ++it1; ++it2;
                            }
                            matrix tmp = matMult(it1_min->m, it2_min->m);
                            if(it1_min->is_tmp) matDestroy(it1_min->m);
                            if(it2_min->is_tmp) matDestroy(it2_min->m);

                            rl.insert(it1_min, data_type{tmp, true});
                            rl.erase(it1_min);
                            rl.erase(it2_min);
                        }
                        res = std::move(rl);
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
                        it_type it1, it2;
                        it1 = it2 = rl.begin();
                        ++it2;
                        matrix tmp = matSum(it1->m, it2->m);
                        if(it1->is_tmp) matDestroy(it1->m);
                        if(it2->is_tmp) matDestroy(it2->m);
                        while(++it2 != rl.end()){
                            matrix aux = matSum(tmp, it2->m);
                            if(it2->is_tmp) matDestroy(it2->m);
                            matDestroy(tmp);
                            tmp = aux;
                        }
                        res.insert(res.begin(), data_type{tmp, true});
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, STAR, ll);
                    matrix tmp = matClos(ll.front().m, 0);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                   // std::cout << "STAR : " << ll.front() << std::endl;
                    res.insert(res.begin(), data_type{tmp, true});
                    break;
                }
                case PLUS:
                {
                    list_type ll;
                    traversal(rpqTree, node->e1, PLUS, ll);
                    matrix tmp = matClos(ll.front().m, 1);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true});
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
                    res.insert(res.begin(), data_type{tmp, true});
                    break;
                }
            }
        }

        void traversal_col_fixed(RpqTree* rpqTree, Tree* node, int parentType, int col, list_type &res){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    if(pred > SIZE){
                        std::cout << "Neg." << std::endl;
                    }else {
                        res.insert(res.begin(), data_type{m_matrices[pred], false});
                    }
                    break;
                }
                case CONC:
                {
                    list_type ll, rl;
                    traversal(rpqTree, node->e1, CONC, ll);
                    traversal_col_fixed(rpqTree, node->e2, CONC, col, rl);
                    rl.splice(rl.begin(), ll);
                    if (parentType == CONC){
                        res = std::move(rl);
                    }else{
                        auto it2 = rl.rbegin(); //last element
                        auto it1 = rl.rbegin(); //last element
                        ++it2; //last element -1
                        matrix tmp = matMult1(it2->m, fullSide,it1->m, col);
                        while(++it2 != rl.rend()){
                            matrix aux = matMult1(it2->m, fullSide,tmp, col);
                            if(it2->is_tmp) matDestroy(it2->m);
                            matDestroy(tmp);
                            tmp = aux;
                        }
                        res.insert(res.begin(), data_type{tmp, true});
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
                        it_type it1, it2;
                        it1 = it2 = rl.begin();
                        ++it2;
                        matrix tmp = matSum1(it1->m, it2->m, fullSide, col);
                        if(it1->is_tmp) matDestroy(it1->m);
                        if(it2->is_tmp) matDestroy(it2->m);
                        while(++it2 != rl.end()){
                            matrix aux = matSum1(tmp, it2->m, fullSide, col);
                            if(it2->is_tmp) matDestroy(it2->m);
                            matDestroy(tmp);
                            tmp = aux;
                        }
                        res.insert(res.begin(), data_type{tmp, true});
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal_col_fixed(rpqTree, node->e1, STAR, col, ll);
                    matrix tmp = matClos1(ll.front().m, 0, fullSide, col);
                    // std::cout << "STAR : " << ll.front() << std::endl;
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true});
                    break;
                }
                case PLUS:
                {
                    list_type ll;
                    traversal_col_fixed(rpqTree, node->e1, PLUS, col, ll);
                    matrix tmp = matClos1(ll.front().m, 1, fullSide, col);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true});
                    break;
                }
                case QUESTION:
                {
                    list_type ll;
                    traversal_col_fixed(rpqTree, node->e1, QUESTION, col, ll);
                    matrix Id = matId(std::max(ll.front().m->height,ll.front().m->width));
                    matrix tmp = matSum1(ll.front().m, Id, fullSide, col);
                    matDestroy(Id);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true});
                    break;
                }
            }
        }

        void traversal_row_fixed(RpqTree* rpqTree, Tree* node, int parentType, int row, list_type &res){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    if(pred > SIZE){
                        std::cout << "Neg." << std::endl;
                    }else {
                        res.insert(res.begin(), data_type{m_matrices[pred], false});
                    }
                    break;
                }
                case CONC:
                {
                    list_type ll, rl;
                    traversal_row_fixed(rpqTree, node->e1, CONC, row, ll);
                    traversal_row_fixed(rpqTree, node->e2, CONC, row, rl);
                    rl.splice(rl.begin(), ll);
                    if (parentType == CONC){
                        res = std::move(rl);
                    }else{
                        rev_it_type it1, it2;
                        it1 = it2 = rl.rbegin(); //last_element
                        ++it2; //last element -1
                        matrix tmp = matMult1(it2->m, row,it1->m, fullSide);
                        while(++it2 != rl.rend()){
                            matrix aux = matMult1(it2->m, row,tmp, fullSide);
                            if(it2->is_tmp) matDestroy(it2->m);
                            matDestroy(tmp);
                            tmp = aux;
                        }
                        res.insert(res.begin(), data_type{tmp, true});
                    }
                    break;
                }

                case OOR:
                {
                    list_type ll, rl;
                    traversal_row_fixed(rpqTree, node->e1, OOR, row, ll);
                    traversal(rpqTree, node->e2, OOR, rl);
                    rl.splice(rl.begin(), ll);
                    if (parentType == OOR){
                        res = std::move(rl);
                    }else{
                        it_type it1, it2;
                        it1 = it2 = rl.begin();
                        ++it2;
                        matrix tmp = matSum1(it1->m, it2->m, row, fullSide);
                        if(it1->is_tmp) matDestroy(it1->m);
                        if(it2->is_tmp) matDestroy(it2->m);
                        while(++it2 != rl.end()){
                            matrix aux = matSum1(tmp, it2->m, row, fullSide);
                            if(it2->is_tmp) matDestroy(it2->m);
                            matDestroy(tmp);
                            tmp = aux;
                        }
                        res.insert(res.begin(), data_type{tmp, true});
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal_row_fixed(rpqTree, node->e1, STAR, row, ll);
                    matrix tmp = matClos1(ll.front().m, 0, row, fullSide);
                    // std::cout << "STAR : " << ll.front() << std::endl;
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true});
                    break;
                }
                case PLUS:
                {
                    list_type ll;
                    traversal_row_fixed(rpqTree, node->e1, PLUS, row, ll);
                    matrix tmp = matClos1(ll.front().m, 1, row, fullSide);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true});
                    break;
                }
                case QUESTION:
                {
                    list_type ll;
                    traversal_row_fixed(rpqTree, node->e1, QUESTION, row, ll);
                    matrix Id = matId(std::max(ll.front().m->height,ll.front().m->width));
                    matrix tmp = matSum1(ll.front().m, Id, row, fullSide);
                    matDestroy(Id);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true});
                    break;
                }
            }
        }

        void traversal_row_col_fixed(RpqTree* rpqTree, Tree* node, int parentType, int row, int col, list_type &res){
            switch(node->type) {
                case STR:{
                    auto pred = rpqTree->getPred(node->pos);
                    if(pred > SIZE){
                        std::cout << "Neg." << std::endl;
                    }else {
                        res.insert(res.begin(), data_type{m_matrices[pred], false});
                    }
                    break;
                }
                case CONC:
                {
                    list_type ll, rl;
                    traversal_row_fixed(rpqTree, node->e1, CONC, row, ll);
                    traversal_col_fixed(rpqTree, node->e2, CONC, col, rl);
                    rl.splice(rl.begin(), ll);
                    if (parentType == CONC){
                        res = std::move(rl);
                    }else{
                        it_type it1, it2;
                        it1 = it2 = rl.begin();
                        ++it2;
                        matrix tmp = matMult1(it2->m, row,it1->m, col);
                        while(++it2 != rl.end()){
                            matrix aux = matMult1(it2->m, row,tmp, col);
                            if(it2->is_tmp) matDestroy(it2->m);
                            matDestroy(tmp);
                            tmp = aux;
                        }
                        res.insert(res.begin(), data_type{tmp, true});
                    }
                    break;
                }

                case OOR:
                {
                    list_type ll, rl;
                    traversal_row_fixed(rpqTree, node->e1, OOR, row, ll);
                    traversal_col_fixed(rpqTree, node->e2, OOR, col, rl);
                    rl.splice(rl.begin(), ll);
                    if (parentType == OOR){
                        res = std::move(rl);
                    }else{
                        it_type it1, it2;
                        it1 = it2 = rl.begin();
                        ++it2;
                        matrix tmp = matSum1(it1->m, it2->m, row, col);
                        if(it1->is_tmp) matDestroy(it1->m);
                        if(it2->is_tmp) matDestroy(it2->m);
                        while(++it2 != rl.end()){
                            matrix aux = matSum1(tmp, it2->m, row, col);
                            if(it2->is_tmp) matDestroy(it2->m);
                            matDestroy(tmp);
                            tmp = aux;
                        }
                        res.insert(res.begin(), data_type{tmp, true});
                    }
                    break;
                }
                case STAR:
                {
                    list_type ll;
                    traversal_row_col_fixed(rpqTree, node->e1, STAR, row, col, ll);
                    matrix tmp = matClos1(ll.front().m, 0, row, col);
                    // std::cout << "STAR : " << ll.front() << std::endl;
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true});
                    break;
                }
                case PLUS:
                {
                    list_type ll;
                    traversal_row_col_fixed(rpqTree, node->e1, PLUS, row, col, ll);
                    matrix tmp = matClos1(ll.front().m, 1, row, col);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true});
                    break;
                }
                case QUESTION:
                {
                    list_type ll;
                    traversal_row_col_fixed(rpqTree, node->e1, QUESTION, row, col, ll);
                    matrix Id = matId(std::max(ll.front().m->height,ll.front().m->width));
                    matrix tmp = matSum1(ll.front().m, Id, row, col);
                    matDestroy(Id);
                    if(ll.front().is_tmp) matDestroy(ll.front().m);
                    res.insert(res.begin(), data_type{tmp, true});
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

        }

        matrix solve_var_to_var(std::string &query){
            list_type res;
            RpqTree rpqTree(query, map_P, SIZE);
            traversal(&rpqTree, rpqTree.root(), STR, res);
            return res.front().m;
        }

        matrix solve_con_to_var(std::string &query, int s_id){
            list_type res;
            RpqTree rpqTree(query, map_P, SIZE);
            traversal_row_fixed(&rpqTree, rpqTree.root(), STR, s_id, res);
            return res.front().m;
        }

        matrix solve_var_to_con(std::string &query, int o_id){
            list_type res;
            RpqTree rpqTree(query, map_P, SIZE);
            traversal_col_fixed(&rpqTree, rpqTree.root(), STR, o_id, res);
            return res.front().m;
        }

        matrix solve_con_to_con(std::string &query, int s_id, int o_id){
            list_type res;
            RpqTree rpqTree(query, map_P, SIZE);
            traversal_row_col_fixed(&rpqTree, rpqTree.root(), STR, s_id, o_id, res);
            return res.front().m;
        }

    };


}

#endif //MATRIX_RPQ_RPQ_SOLVER_HPP