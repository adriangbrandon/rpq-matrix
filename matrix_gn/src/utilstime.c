//
// Created by Adrián on 5/5/23.
//

#include "utilstime.h"

uint64_t user_now(){
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    return (r_usage.ru_utime.tv_sec *1000000 + r_usage.ru_utime.tv_usec)*1000;
}

uint64_t system_now(){
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    return (r_usage.ru_stime.tv_sec * 1000000 + r_usage.ru_stime.tv_usec)*1000;
}

uint64_t elapsed_now(){
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    return ((r_usage.ru_stime.tv_sec + r_usage.ru_utime.tv_sec) * 1000000
            + r_usage.ru_utime.tv_usec + r_usage.ru_stime.tv_usec)*1000;
}
