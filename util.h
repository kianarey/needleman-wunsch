#ifndef UTIL_H__
#define UTIL_H__

#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#define TIME_IT(ROUTINE_NAME__, LOOPS__, ACTION__)\
{\
    printf("    Timing '%s' started\n", ROUTINE_NAME__);\
    struct timeval tv;\
    struct timezone tz;\
    const clock_t startTime = clock();\
    gettimeofday(&tv, &tz); long GTODStartTime =  tv.tv_sec * 1000 + tv.tv_usec / 1000 ;\
    for (int loops = 0; loops < (LOOPS__); ++loops)\
    {\
        ACTION__;\
    }\
    gettimeofday(&tv, &tz); long GTODEndTime =  tv.tv_sec * 1000 + tv.tv_usec / 1000 ;\
    const clock_t endTime = clock();\
    const clock_t elapsedTime = endTime - startTime;\
    const double timeInSeconds = (elapsedTime/(double)CLOCKS_PER_SEC);\
    printf("        GetTimeOfDay Time (for %d iterations) = %g\n", LOOPS__, (double)(GTODEndTime - GTODStartTime) / 1000. );\
    printf("        Clock Time        (for %d iterations) = %g\n", LOOPS__, timeInSeconds );\
    printf("    Timing '%s' ended\n", ROUTINE_NAME__);\
}

void print_score_matrix(int *score_matrix, int rows, int cols);

#endif