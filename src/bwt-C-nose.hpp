/*
 * bwt_nose.hpp
 * Copyright (C) 2020 Author removed for double-blind evaluation
 * 
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


#ifndef BWT_NOSE_T
#define BWT_NOSE_T

#include "Config.hpp"
#include "wt_intersection.hpp"

using namespace std;


class bwt_nose
 {
    bwt_type L;

    C_type *C_bv;
    C_rank_type C_rank;
    C_select_type C_select;
    C_select0_type C_select0;
    std::vector<uint64_t> C_a;

    public:

        //std::vector<uint64_t> C;
		
        bwt_nose() {;}

	bwt_nose(int_vector<> & _L, vector<uint64_t> & _C)  
        {
            construct_im(L, _L);
            //C_a = _C;
            //bit_vector C_aux = bit_vector(_C[_C.size()-1]+1+_C.size(), 0);

            for (uint64_t i=0; i < _C.size(); i++)
                C_a.push_back(_C[i]);
         

            //C_bv = new C_type(C_aux);
            //util::init_support(C_rank,C_bv);
            //util::init_support(C_select,C_bv);
            //util::init_support(C_select0,C_bv);
        }

        ~bwt_nose() {;}

        uint64_t operator[](uint64_t i) 
	{
	    return L[i];
	}


        uint64_t size()
        {
            return sdsl::size_in_bytes(L) + sizeof(uint64_t)*C_a.size(); //sdsl::size_in_bytes(*C_bv) 
                 //+ sdsl::size_in_bytes(C_rank) + sdsl::size_in_bytes(C_select)
                 //+ sdsl::size_in_bytes(C_select0);
        }

	void save(string filename)
        {
            sdsl::store_to_file(L, filename+".L");

            //sdsl::store_to_file(*C_bv, filename+".C");
            std::ofstream fp_C(filename+".Ca");
            for (uint64_t i = 0; i < C_a.size(); i++) {
                fp_C << C_a[i] << endl;
            }
            
        }

        void load(string filename)
        {
//            cout << "OK, llego al load" << endl; fflush(stdout);
            sdsl::load_from_file(L, filename+".L");

//            cout << "OK, ya cargo el L_P" << endl; fflush(stdout);
            std::ifstream fp_C(filename+".Ca");
/*            C_bv = new C_type; //bit_vector;
            sdsl::load_from_file(*C_bv, filename+".C");
            util::init_support(C_rank,C_bv);
            util::init_support(C_select,C_bv);
            util::init_support(C_select0,C_bv);
*/
/*            cout << "Ahora calcula el rank" << endl; fflush(stdout);
            uint64_t r = C_rank(C_bv->size()-1), t;

            cout << r << " ones in C_bv" << endl; fflush(stdout);

            for (uint64_t i = 0; i < C_bv->size(); i++ ) {
                //cout << i << endl; fflush(stdout);
                if ((*C_bv)[i]) {*/
//                   C_a.push_back(i/*C_select(i)*/);
//                   cout << C_a.size() << endl;
//                }
//            }
              
/*
            for (uint64_t i = 0; i < C_a.size(); i++) {
                cout << i << " " << C_a[i] << endl; fflush(stdout);
                fp_C << C_a[i] << endl;
            }
            exit(1);*/
            uint64_t t, i=0;

            while (!fp_C.eof()) {
                fp_C >> t;
                //t -= i; i++;
                C_a.push_back(t);
            }

            /*for (i = 0; i < C_a.size()-1; i++) {
	        get_C(i);
                //cin >> t;
	    }*/
            //exit(1);



        }

        uint64_t LF(uint64_t i)
        {
            uint64_t s = L[i];		 	
            return get_C(s) + L.rank(i, s) - 1;
        }

        uint64_t nElems(uint64_t val)
        {
            return get_C(val+1)-get_C(val);
        }

        pair<uint64_t, uint64_t> 
        backward_step(uint64_t left_end, uint64_t right_end, uint64_t value)
        {
            uint64_t s = L.rank(left_end, value);
            uint64_t e = L.rank(right_end+1, value)-1;
            return pair<uint64_t, uint64_t>(s, e);
        }

        pair<uint64_t, uint64_t>
        backward_step_test(uint64_t left_end, uint64_t right_end, uint64_t value, 
                           uint64_t s, uint64_t e)
        {
            return pair<uint64_t, uint64_t>(s, e-1);
        }

        /*inline uint64_t bsearch_C(uint64_t value)
        {
            uint64_t r = C_rank(C_select0(value+1));
            return r;  
        }*/

        inline uint64_t get_C(uint64_t v) const
        {   
            //if (C_select(v+1)-v != C_a[v])
            //    cout << "Deberia ser = " << C_select(v+1)-v << "  pero es = " << C_a[v] << endl;
            return C_a[v]; // OJO aqui, posible punto de error
        }

        inline uint64_t ranky(uint64_t pos, uint64_t val)
        {
            return L.rank(pos, val);
        }

        inline uint64_t rank(uint64_t pos, uint64_t val)
        {
            return L.rank(get_C(pos), val);
        }

        inline uint64_t select(uint64_t _rank, uint64_t val)
        {
            return L.select(_rank, val);
        }

        inline std::pair<uint64_t,uint64_t> select_next(uint64_t pos, uint64_t val, uint64_t n_elems)
        {
            return L.select_next(get_C(pos), val, n_elems);
        }

        inline uint64_t min_in_range(uint64_t l, uint64_t r)
        {
            return L.range_minimum_query(l, r);
        }

        inline uint64_t range_next_value(uint64_t x, uint64_t l, uint64_t r)
        {
            return L.range_next_value(x, l, r);
        }

        inline std::pair<uint64_t, uint64_t> range_next_value_pos(uint64_t x, uint64_t l, uint64_t r)
        {
            return L.range_next_value_pos(x, l, r);
        }
        std::vector<pair<uint64_t, uint64_t>>
        //inline void
        values_in_range(uint64_t pos_min, uint64_t pos_max, uint64_t sigma/*, std::vector<uint64_t> & values, uint64_t & k*/)
        { 
            //interval_symbols(L, pos_min, pos_max+1, k, values, r_i, r_j);
            return L.range_search_2d(pos_min, pos_max, 1, sigma).second;
        } 

        // backward search for pattern of length 1
        pair<uint64_t, uint64_t> backward_search_1_interval(uint64_t P) const
        {
            return pair<uint64_t, uint64_t>(get_C(P), get_C(P+1)-1);
        }

        // backward search for pattern of length 1
        pair<uint64_t, uint64_t> backward_search_1_rank(uint64_t P, uint64_t S) const
        {
            uint64_t s = L.rank(get_C(P), S);
            uint64_t e = L.rank(get_C(P+1), S);
            return pair<uint64_t, uint64_t>(s, e);
        }

	// backward search for pattern PQ of length 2
	// returns an empty interval if search is unsuccessful
        pair<uint64_t, uint64_t> 
        backward_search_2_interval(uint64_t P, pair<uint64_t, uint64_t> &I) const
        {
            return pair<uint64_t, uint64_t>(get_C(P)+I.first, get_C(P)+I.second-1);
        }

        pair<uint64_t, uint64_t> 
        backward_search_2_rank(uint64_t P, uint64_t S, pair<uint64_t, uint64_t> &I) const
        {
            uint64_t c = get_C(P);
            uint64_t s = L.rank(c+I.first, S);
            uint64_t e = L.rank(c+I.second, S);
            return pair<uint64_t, uint64_t>(s, e);			
        }

        std::vector<uint64_t>
        all_values_in_range(uint64_t lb, uint64_t rb) const
        {
            return L.all_values_in_range(lb, rb);
        }




        template<typename word_t>
        void mark(uint64_t c, std::vector<word_t>& B_wt, word_t B_c) const
        {
	    L.mark<word_t>(c, B_wt, B_c);
	
	} 

        template<typename word_t>
        void unmark(uint64_t c, std::vector<word_t>& B_wt) const
        {
            L.unmark<word_t>(c, B_wt);

        }


        template<typename word_t>
        std::vector<uint64_t> 
        all_active_p_values_in_range(uint64_t lb, uint64_t rb, 
                                     std::vector<word_t>& B_wt, word_t D,
                                     bool report=true) const 
        {
            return L.all_active_p_values_in_range<word_t>(lb, rb, B_wt, D);
        }

        template<typename word_t>
        //std::vector<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>>
        void
        all_active_p_values_in_range_test(uint64_t lb, uint64_t rb,
                                     std::vector<word_t>& B_wt, word_t D,
                                     std::vector<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>>& res_vec,
                                     bool report=true) const
        {
            /*return*/ L.all_active_p_values_in_range_test<word_t>(lb, rb, B_wt, D, res_vec);
        }

        template <typename word_t>
        std::vector<std::pair<uint64_t, word_t>>
        all_active_s_values_in_range(uint64_t lb, uint64_t rb,
                                   initializable_array<word_t>& D_wt, word_t D,
                                   bool report=true) const
        {
	    return L.all_active_s_values_in_range(lb, rb, D_wt, D); 	
	} 

        template <typename word_t>
        //std::vector<std::tuple<uint64_t, word_t,std::pair<uint64_t,uint64_t>>>
        void
        all_active_s_values_in_range_test(uint64_t lb, uint64_t rb,
                                   initializable_array<word_t>& D_wt, word_t D,
                                   std::vector<std::tuple<uint64_t, word_t,std::pair<uint64_t,uint64_t>>>& res_vec,
                                   bool report=true) const
        {
            /*return*/ L.all_active_s_values_in_range_test(lb, rb, D_wt, D, res_vec);
        }


        void all_values_in_range_bounded(uint64_t lb, uint64_t rb, 
                                         std::vector<uint64_t>& res_vec, uint64_t bound, 
     		         	         bool report=true) const 
        {
	    L.all_values_in_range_bounded(lb, rb, res_vec, bound, report);
	};

        
	std::vector<std::pair<uint64_t,uint64_t>>
	intersect(const std::vector<std::array<uint64_t, 2ul>>& ranges)
	{
            return sdsl::intersect<bwt_type>(L, ranges);
        };

    std::vector<uint64_t>
    intersect_nofreq(const std::vector<std::array<uint64_t, 2ul>>& ranges)
    {
        return sdsl::intersect_nofreq<bwt_type>(L, ranges);
    };

 };
#endif
