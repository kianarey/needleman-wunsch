#ifndef NEEDLEMAN_WUNSCH_CPU_H__
#define NEEDLEMAN_WUNSCH_CPU_H__

#define MISMATCH -3;
#define MATCH 1;
#define GAP -2;

using namespace std;
#include <string>

int NeedleMan_Munsch(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols);

void traceBack(int *scoreMatrix, int rows, int cols, char *seq1, char* seq2);


#endif