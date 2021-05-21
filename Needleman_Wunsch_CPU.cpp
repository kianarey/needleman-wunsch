#include "Needleman_Wunsch_CPU.h"
#include "util.h"
#include <string>
#include <fstream>
#include <stdio.h>
#include <string.h>


using namespace std;

/// <summary>
/// Function that performs needleman wunsch alg to align 2 sequences,
/// updates given score matrix and returns final alignment score
/// </summary>
/// <param name="scoreMatrix">pointer to flattened 2D array to save score calculations to</param>
/// <param name="seq1">first nucleotide sequence to align</param>
/// <param name="seq2">second nucleotide sequence to align</param>
/// <returns>final alignment score of 2 sequences</returns>
int NeedleMan_Munsch(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols) {
	// rows correspond to seq1, cols to seq2
	// set gap alignment values for row 0
	for (int x = 0; x < cols; x++) {
		scoreMatrix[x] = x * GAP;
	}

	// set gap alignment values for col 0
	for (int x = 0; x < rows * cols; x += cols) {
		scoreMatrix[x] = (x / cols) * GAP;
	}

	// fill in matrix row by row starting at row 1
	for (int curRow = 1; curRow < rows; curRow++) {
		for (int curCol = 1; curCol < cols; curCol++) {
			// calculate score from going down, this introduces a gap
			int	down = scoreMatrix[(curRow - 1) * cols + curCol] + GAP;

			// calculate score from going right, this introduces a gap
			int	right = scoreMatrix[curRow * cols + (curCol - 1)] + GAP;

			// calculate score from going diagonal
			int	diag = 0;
			if (seq1[curRow - 1] == seq2[curCol - 1]) {
				//seq1 row and seq2 col nucleotides match
				diag = scoreMatrix[(curRow - 1) * cols + (curCol - 1)] + MATCH;
			}
			else {
				//seq1 row and seq2 col nucleotides do not match
				diag = scoreMatrix[(curRow - 1) * cols + (curCol - 1)] + MISMATCH;
			}

			// update cell value scoreand trace values based on highest score
			if (diag >= down && diag >= right) {
				//diag produced best score
				scoreMatrix[curRow * cols + curCol] = diag;
			}
			else if (right >= diag && right >= down) {
				// right produced best score
				scoreMatrix[curRow * cols + curCol] = right;
			}
			else {
				// down produced best score
				scoreMatrix[curRow * cols + curCol] = down;
			}
		}
	}

	// return final score
	return scoreMatrix[(rows - 1) * cols + (cols - 1)];
}



void traceBack(int *scoreMatrix, int rows, int cols, char *seq1, char* seq2) {
    string seq1Aligned = "";
    string seq2Aligned = "";

    // Start from last cell of the table
    int curRow = rows - 1;
    int curCol = cols - 1;

    // Run the traceback until it reaches cell (0,0)
	// Rows correspond to seq1, cols to seq2
    while (curRow > 0 || curCol > 0) {
        if(curRow == 0) {
            // At the end of seq 1
            seq1Aligned = "-" + seq1Aligned;
            seq2Aligned = seq2[curCol] + seq2Aligned;
            curCol--;
        } else if (curCol == 0) {
            // At the end of seq 2
            seq1Aligned = seq1[curRow] + seq1Aligned;
            seq2Aligned = "-" + seq2Aligned;
            curRow--;
        } else {
            // Find the reverse direction
            int	up = scoreMatrix[(curRow - 1) * cols + curCol] + GAP;
            int left = scoreMatrix[(curRow) * cols + curCol - 1] + GAP;

            if (scoreMatrix[curRow*cols + curCol] == up) {
                // If prev cell = up :
                //  add gap to seq 2, copy seq 1's character
                seq1Aligned = seq1[curRow] + seq1Aligned;
                seq2Aligned = "-" + seq2Aligned;
                curRow--;
            } else if (scoreMatrix[curRow*cols + curCol] == left) {
                // If prev cell = left :
                //  add gap to seq 1, copy seq 2's character
                seq1Aligned = "-" + seq1Aligned;
                seq2Aligned = seq2[curCol] + seq2Aligned;
                curCol--;
            } else {
                // prev cell = diag
                //  copy both characters, move to prev cell
                seq1Aligned = seq1[curRow] + seq1Aligned;
                seq2Aligned = seq2[curCol] + seq2Aligned;
                curRow--;
                curCol--;
            }
        }
        
    }
    
    // Print aligned seq 1 to output file.
    ofstream outputSeq1("Seq1Aligned.txt");
    char * seq1Aligned_ = (char *)seq1Aligned.c_str();
    outputSeq1.write(seq1Aligned_, strlen(seq1Aligned_));
    outputSeq1.close();

    // print aligned seq 2 to output file
    ofstream outputSeq2("Seq2Aligned.txt");
    char * seq2Aligned_ = (char *)seq2Aligned.c_str();
    outputSeq2.write(seq2Aligned_, strlen(seq2Aligned_));
    outputSeq2.close();

}
