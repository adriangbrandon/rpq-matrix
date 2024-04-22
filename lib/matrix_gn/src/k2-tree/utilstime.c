//
// Created by Adrián on 5/5/23.
//

#include "k2-tree/utilstime.h"


void time_beg() {
    clock_gettime(CLOCK_REALTIME, &beg);
}

double time_diff(){
    struct timespec end;
    clock_gettime(CLOCK_REALTIME, &end);
    double start_sec, end_sec, elapsed_sec;
    start_sec = beg.tv_sec + beg.tv_nsec / NANO_PER_SEC;
    end_sec = end.tv_sec + end.tv_nsec / NANO_PER_SEC;
    elapsed_sec = end_sec - start_sec;
    return elapsed_sec;
}
