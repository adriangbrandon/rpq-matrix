//
// Created by Adrián on 16/5/23.
//

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

std::string file_name(const uint i, const std::string &index){
    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << i;
    return index + "/" + ss.str() + ".mat";
}