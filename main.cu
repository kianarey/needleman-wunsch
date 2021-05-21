#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <time.h> 
#include "util.h"

#include "Needleman_Wunsch_CPU.h"
#include "NW_Single_Kernel.cu"

using namespace std;

// Function declarations
char RandSeqLetter();
void GenerateRandomInputSequence(char * sequence, int size);


int main(int argc, char* argv[])
{
    char * inputSequence1;
    char * inputSequence2;
    int seq1Size;
    int seq2Size;
    srand (time(NULL));
    int n = 10; // number of iterations
    int kernel_mode = 0; // 1 = launch_diag_block_kernel, 2 = launch_multi_kernel_diag, 3 = launch_single_block_kernel

    /* 
        The command line for calling your program should be one of the following forms: 
        programname [seq1.fna] [seq2.fna] [kernel_mode]
        programname [sequence size] [kernel_mode]
        programname
    */
    if (argc == 4){
        // Two input files for sequences to align are specified
        string inputText1;
        string inputText2;

        ifstream inputFile1(argv[1]);
        getline (inputFile1, inputText1);
        inputSequence1  = (char *) inputText1.c_str();
        
        ifstream inputFile2(argv[2]);
        getline (inputFile2, inputText2);
        inputSequence2  = (char *) inputText2.c_str();

        kernel_mode = atoi(argv[3]);

    }
    else if (argc == 3) {
        // Two sequences of specified size will be aligned
        seq1Size = atoi(argv[1]);
        seq2Size = atoi(argv[1]);
        inputSequence1 = (char*)malloc(sizeof(char) * seq1Size);
        inputSequence2 = (char*)malloc(sizeof(char) * seq2Size);

        GenerateRandomInputSequence(inputSequence1, seq1Size);
        GenerateRandomInputSequence(inputSequence2, seq2Size);

        kernel_mode = atoi(argv[2]);
    }
    else {
        // Two sequences of random size (1-100) will be aligned
        seq1Size = rand() % 100 + 1;
        seq2Size = rand() % 100 + 1;
        
        inputSequence1 = (char*) malloc(sizeof(char) * seq1Size);
        inputSequence2 = (char*) malloc(sizeof(char) * seq2Size);
        
        GenerateRandomInputSequence(inputSequence1, seq1Size);
        GenerateRandomInputSequence(inputSequence2, seq2Size);

        kernel_mode = 1; // default if not specified
    }

    int rows = seq1Size + 1;
    int cols = seq2Size + 1;

    int* scoreMatrix = (int *)malloc(sizeof(int)*rows*cols);
    int* h_scoreMatrix = (int*)malloc(sizeof(int) * (rows * cols));

    // Perform CPU calculation
    TIME_IT("Needleman_Wunsch_CPU", n, NeedleMan_Munsch(scoreMatrix, inputSequence1, inputSequence2, rows, cols);)
    
    // Perform GPU calculation
    computeOnDevice(h_scoreMatrix, inputSequence1, inputSequence2, rows, cols, n, kernel_mode);

    // Compare score matrix of CPU and GPU
    int passed = 1;
    for (int i = 0; i < (rows * cols); i++) {
        if (h_scoreMatrix[i] != scoreMatrix[i]) {
            printf("GPU score[%i]: %i, CPU score[%i]: %i\n", i, h_scoreMatrix[i], i, scoreMatrix[i]);
            print_score_matrix(h_scoreMatrix, rows, cols);
            passed = 0;
            break;
        }
    }
    //  printf("Rows: %i, Cols: %i\n", rows, cols);
    (passed) ? printf("\n    Test PASSED\n") : printf("\n    Test FAILED\n");

    // Traces back the scorematrix
    // Creates aligned sequence and prints them to output files
    traceBack(h_scoreMatrix, rows, cols, inputSequence1, inputSequence2);

    // Clean up memory
    free(scoreMatrix);
    free(h_scoreMatrix);
    scoreMatrix = NULL;
    return 0;
}

void GenerateRandomInputSequence(char * sequence, int size) {
    for(int i=0; i< size; i++) {
        sequence[i] = RandSeqLetter();
    }
}

char RandSeqLetter() {
    int num = rand() % 4;
    switch (num) {
    case 0:
        return 'A';
    case 1:
        return 'G';
    case 2:
        return 'C';
    case 3:
        return 'T';
    }

    //should never reach here
    return 'A';
}

