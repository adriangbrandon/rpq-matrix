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
// Created by Adrián on 31/12/20.
//

#ifndef UTILITIES_FILE_HPP
#define UTILITIES_FILE_HPP


#include <vector>
#include <iostream>
#include <fstream>
#include <cctype>
#include <sys/stat.h>
#include <dirent.h>
#include <stdint.h>


namespace util {


    namespace file {

        std::vector<std::string> read_directory(const std::string& name)
        {
            struct dirent **namelist;
            int n;
            std::vector<std::string> v;
            n = scandir(name.c_str(), &namelist, 0, alphasort);
            if (n > 0){
                while (n--) {
                    if(namelist[n]->d_type == DT_REG){
                        v.push_back(namelist[n]->d_name);
                    }
                    free(namelist[n]);
                }
                free(namelist);
            }
            return v;
        }

        bool create_directory(const std::string& name){
            int status;
            status = mkdir(name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            return (status == 0);
        }

        uint64_t file_size(const std::string& file){
            struct stat fs;
            if(stat(file.c_str(), &fs) != 0) {
                return 0;
            };
            return fs.st_size;
        }

        bool file_exists(const std::string &file){
            struct stat fs;
            return (stat(file.c_str(), &fs) == 0);
        }

        bool remove_file(const std::string &file){
            return (remove(file.c_str())==0);
        }

        std::string remove_path(const std::string &file){
            auto pos_last_slash = file.find_last_of("\\/");
            if(pos_last_slash != -1) {
                return file.substr(pos_last_slash+1);
            }
            return file;
        }

        bool end_slash(const std::string &file){
            auto pos_last_slash = file.find_last_of("\\/");
            return (pos_last_slash == file.length()-1);
        }

        std::string remove_extension(const std::string &file){
            auto pos_last_slash = file.find_last_of("\\/");
            auto pos_last_point = file.find_last_of('.');
            if(pos_last_point != -1 && (pos_last_slash == -1 || pos_last_slash < pos_last_point)){
                return file.substr(0, pos_last_point);
            }
            return file;
        }

        std::string index_file(const std::string &index_name, const char* argv[], const size_t length){
            std::string result = index_name;
            if(length > 1){
                std::string dataset_name = remove_extension(remove_path(argv[1]));
                result += "_" + dataset_name;
            }
            for(size_t i = 2; i < length; ++i){
                result += "_" + std::string(argv[i]);
            }
            return result;
        }


        template<class t_value>
        void write_to_file(const std::string& file, const std::vector<t_value> &container){
            std::ofstream out(file, std::ios::out | std::ios::binary );
            out.write((char*) &container[0], container.size() * sizeof(t_value));
        }

        template<class t_value>
        void read_from_file(const std::string& file, std::vector<t_value> &container){
            std::ifstream in(file, std::ios::in | std::ios::binary );
            uint64_t size = file_size(file);
            container.resize(size / sizeof(t_value));
            in.read((char*) &container[0], container.size() * sizeof(t_value));
        }
    }



}



#endif //UTILITIES_FILE_HPP
