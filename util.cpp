#include "util.h"

void print_score_matrix(int *score_matrix, int rows, int cols){
    int i, j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++)
            printf("%d ",score_matrix[i*cols + j]);
        printf("\n");
    }
}