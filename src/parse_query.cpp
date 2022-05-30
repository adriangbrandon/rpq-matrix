/*
 * triple_bwt_rpq.hpp
 * Copyright (C) 2020 Author name removed for double blind evaluation
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

#ifndef PARSE_QUERY
#define PARSE_QUERY

#include <cstdint>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <stdio.h>
#include <stdlib.h>


std::string parse(const std::string& rpq, int64_t& i,
                  std::unordered_map<std::string, uint64_t>& predicates_map, uint64_t max_P)
{
    std::string str_aux;
    char _operator;

    while (i < rpq.size()) {
        if (rpq[i]==')')
            return str_aux;
        else if (rpq[i] == '<') {
            std::string s("<");
            i++;
            if (rpq[i]=='%') {
                i++;
                while (i < rpq.size() and rpq[i]!='>') {
                    s += rpq[i];
                    i++;
                }
                str_aux += s + '>';
                i++;
            } else {
                std::string s_aux("<");
                s += '%';
                while (i < rpq.size() and rpq[i]!='>') {
                    s += rpq[i];
                    s_aux += rpq[i];
                    i++;
                }
                s += '>';
                s_aux += '>';
                predicates_map[s] = max_P + predicates_map[s_aux];
                str_aux += s;
                i++;
            }
        } else if (rpq[i]=='(') {
                i++;
                std::string s = parse(rpq, i, predicates_map, max_P);
                //std::reverse(s.begin(), s.end());
                str_aux += '(' + s + ')';
                i++;
        } else {
                str_aux += rpq[i];
                i++;
        }
    }
    return str_aux;
};


std::string parse_reverse(const std::string& rpq, int64_t& i, 
                          std::unordered_map<std::string, uint64_t>& predicates_map, uint64_t max_P)
{
    std::string str_aux;
    char _operator;

    while (i >= 0) {
        if (rpq[i]=='(')
	    return str_aux;
	else if (rpq[i]=='*' or rpq[i]=='+' || rpq[i]=='?') {
            _operator = rpq[i];

            if (rpq[i-1]==')') {
                i -= 2;
                std::string s = parse_reverse(rpq, i, predicates_map, max_P);
                //std::reverse(s.begin(), s.end());
                str_aux += '(' + s + ')' + _operator;
                i--;
            } else {
                std::string s(">");
                std::string s_aux(">");
                i--;
                while (i >= 0 and rpq[i]!='<') {
                    s += rpq[i];
                    if (rpq[i] != '%') 
                        s_aux += rpq[i];
                    i--;
                }
                s += '<';
                s_aux += '<';
                std::reverse(s.begin(), s.end());
                if (s[1] == '%') {
                    std::reverse(s_aux.begin(), s_aux.end());
                    predicates_map[s] = max_P + predicates_map[s_aux];
                }
                str_aux += s + _operator;
                i--;
            }
        } else if (rpq[i] == '>') {
            std::string s(">");
            std::string s_aux(">"); 
            i--;
            while (i >= 0 and rpq[i]!='<') {
                s += rpq[i];
                if (rpq[i] != '%')
                    s_aux += rpq[i];
                i--;
            }
            s += '<';
            s_aux += '<';
            std::reverse(s.begin(), s.end());
            if (s[1] == '%') {
                std::reverse(s_aux.begin(), s_aux.end());
                predicates_map[s] = max_P + predicates_map[s_aux];
            }
            str_aux += s;
            i--;
        } else if (rpq[i]==')') {
            --i;
            std::string s = parse_reverse(rpq, i, predicates_map, max_P);
            //std::reverse(s.begin(), s.end());
            str_aux += '(' + s + ')';
            i--;
        } else {
            str_aux += rpq[i];
            i--;
        }
    }
    return str_aux;
};

#endif
