#ifndef _NW_SINGLE_KERNEL_H_
#define _NW_SINGLE_KERNEL_H_

#define BLOCK_SIZE 32

#include "Needleman_Wunsch_CPU.h"

__global__ void nw_single_kernel(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols, int startCol, int startRow) {

    // Calculate row index of scoreMatrix element
    int curRow = blockIdx.y * blockDim.y + threadIdx.y + startRow;

    // Calculate column index of scoreMatrix element
    int curCol = blockIdx.x * blockDim.x + threadIdx.x + startCol;

	int startDiag = startCol + startRow;
	int endDiag = startDiag + 2 * BLOCK_SIZE - 1;

    for (int curDiag = startDiag; curDiag < endDiag; curDiag++) {
        // sequentially calculate each element in diagonals
        __syncthreads();

        // verify element is part of diagonal and calculate score
        if (curRow + curCol == curDiag && curRow < rows && curCol < cols && curRow > 0 && curCol > 0 && curRow < startRow + BLOCK_SIZE && curCol < startCol + BLOCK_SIZE) {
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
}

__global__ void initialize_matrix(int* scoreMatrix, int rows, int cols) {

	int index = blockDim.x * blockIdx.x + threadIdx.x;

	// set gap alignment values for row 0
	if (index < cols) {
		scoreMatrix[index] = index * GAP;
	}

	__syncthreads();

	// set gap alignment values for col 0
	if (index < rows) {
		scoreMatrix[index * cols] = index * GAP;
	}

}

__global__ void nw_single_diag_kernel(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols, int blockDiag, int blockCols) {

	// Value based on top left cell in block
	int startRow;
	int startCol;
	if (blockDiag <= blockCols) {
		startRow = blockIdx.y * blockDim.y;
		startCol = (blockDiag - 1 - blockIdx.y) * blockDim.x;
	}
	else {
		// requires additional offset
		startRow = (blockIdx.y + blockDiag - blockCols) * blockDim.y;
		startCol = (blockCols - 1 - blockIdx.y) * blockDim.x;
	}

	// Calculate row index of scoreMatrix element
	int curRow = threadIdx.y + startRow;

	// Calculate column index of scoreMatrix element
	int curCol = threadIdx.x + startCol;

	int startDiag = startCol + startRow;
	int endDiag = startDiag + blockDim.x + blockDim.y - 1;

	for (int curDiag = startDiag; curDiag < endDiag; curDiag++) {
		// sequentially calculate each element in diagonals
		__syncthreads();

		// verify element is part of diagonal and calculate score
		if (curRow + curCol == curDiag && curRow < rows && curCol < cols && curRow > 0 && curCol > 0 && curRow < startRow + blockDim.y && curCol < startCol + blockDim.y) {
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
}



__global__ void nw_single_diag_kernel_shared(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols, int blockDiag, int blockCols) {

	// every thread in the block will place an element into seq1_shared and seq2_shared
	__shared__ char seq1_shared[BLOCK_SIZE*BLOCK_SIZE];
	__shared__ char seq2_shared[BLOCK_SIZE*BLOCK_SIZE];

	// Value based on top left cell in block
	int startRow;
	int startCol;
	if (blockDiag <= blockCols) {
		startRow = blockIdx.y * blockDim.y;
		startCol = (blockDiag - 1 - blockIdx.y) * blockDim.x;
	}
	else {
		// requires additional offset
		startRow = (blockIdx.y + blockDiag - blockCols) * blockDim.y;
		startCol = (blockCols - 1 - blockIdx.y) * blockDim.x;
	}

	// Calculate row index of scoreMatrix element
	int curRow = threadIdx.y + startRow;

	// Calculate column index of scoreMatrix element
	int curCol = threadIdx.x + startCol;

	int startDiag = startCol + startRow;
	int endDiag = startDiag + 2 * BLOCK_SIZE - 1;

	// Variables that flag whether the current block is within the first row/column of blocks
	// This becomes important when initializing the matrix
	int top_block_row;
	int left_block_col;

	// These are used to compare elements between the sequences in later kernel code
	char seq1_char;
	char seq2_char;

	if(curRow < rows){
		if(startRow != 0){
			seq1_shared[threadIdx.x + threadIdx.y*BLOCK_SIZE] = seq1[threadIdx.x + threadIdx.y*BLOCK_SIZE + startRow - 1];
			top_block_row = 0;
		}
		else{
			seq1_shared[threadIdx.x + threadIdx.y*BLOCK_SIZE] = seq1[threadIdx.x + threadIdx.y*BLOCK_SIZE + startRow];
			top_block_row = 1;
		}
	}
	if(curCol < cols){
		if(startCol != 0){
			seq2_shared[threadIdx.x + threadIdx.y*BLOCK_SIZE] = seq2[threadIdx.x + threadIdx.y*BLOCK_SIZE + startCol - 1];
			left_block_col = 0;
		}
		else{
			seq2_shared[threadIdx.x + threadIdx.y*BLOCK_SIZE] = seq2[threadIdx.x + threadIdx.y*BLOCK_SIZE + startCol];
			left_block_col = 1;
		}
	}

	__syncthreads();


	for (int curDiag = startDiag; curDiag < endDiag; curDiag++) {
		// sequentially calculate each element in diagonals
		__syncthreads();

		// verify element is part of diagonal and calculate score
		if (curRow + curCol == curDiag && curRow < rows && curCol < cols && curRow > 0 && curCol > 0 && curRow < startRow + BLOCK_SIZE && curCol < startCol + BLOCK_SIZE) {
			// calculate score from going down, this introduces a gap
			int	down = scoreMatrix[(curRow - 1) * cols + curCol] + GAP;

			// calculate score from going right, this introduces a gap
			int	right = scoreMatrix[curRow * cols + (curCol - 1)] + GAP;

			// calculate score from going diagonal
			int	diag = 0;

			// if it's not a block on the first row
			if(top_block_row == 0)
				seq1_char = seq1_shared[threadIdx.y];
			else
				seq1_char = seq1_shared[threadIdx.y - 1];

			if(left_block_col == 0)
				seq2_char = seq2_shared[threadIdx.x];
			else
				seq2_char = seq2_shared[threadIdx.x - 1];

			if (seq1_char == seq2_char) {
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
}

void launch_single_block_kernel(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols) {
	// fill in first row & col of the matrix
	initialize_matrix << < ceil(max(cols, rows) / (float)256), 256 >> > (scoreMatrix, rows, cols);

    // requires 1 thread/element
	// 1 kernel launch per block (sequentially calculates blocks)
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE, 1);
	for (int startRow = 1; startRow < rows; startRow += BLOCK_SIZE) {
		for (int startCol = 1; startCol < cols; startCol += BLOCK_SIZE) {
			nw_single_kernel << < 1, dimBlock >> > (scoreMatrix, seq1, seq2, rows, cols, startCol, startRow);
		}
	}
}

void launch_diag_block_kernel(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols) {
	// fill in first row & col of the matrix
	initialize_matrix << < ceil(max(cols, rows) / (float)256), 256 >> > (scoreMatrix, rows, cols);

	// requires 1 thread/element
	// 1 kernel launch per diagonal of blocks (sequentially calculates diagonal)
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE, 1);
	int blockRows = ceil((float)rows / BLOCK_SIZE);
	int blockCols = ceil((float)cols / BLOCK_SIZE);
	for (int diag = 1; diag < blockRows + blockCols ; diag++) {
		// determine number of blocks along diagonal
		int diagBlocks;
		if (blockRows <= blockCols) {
			if (diag <= blockRows) {
				diagBlocks = diag;
			}
			else if (diag <= blockCols) {
				diagBlocks = blockRows;
			}
			else {
				diagBlocks = blockRows - (diag - blockCols);
			}
		}
		else {
			if (diag <= blockCols) {
				diagBlocks = diag;
			}
			else if (diag <= blockCols) {
				diagBlocks = blockCols;
			}
			else {
				diagBlocks = blockCols - (diag - blockRows);
			}
		}
		dim3 dimGrid(1, diagBlocks, 1);
		nw_single_diag_kernel_shared << < dimGrid, dimBlock >> > (scoreMatrix, seq1, seq2, rows, cols, diag, blockCols);
	}
}

/*
======================================
	Multi Kernel Implementation / Start
======================================
*/


// Utility function to compute a given cell and write result to scoreMatrix[cur_row *cols + cur_col]
__device__ void compute_cell(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols, int cur_row, int cur_col ) {


	// // calculate score from going down, this introduces a gap
	int	down = scoreMatrix[(cur_row - 1) * cols + cur_col] + GAP;

	// // calculate score from going right, this introduces a gap
	int	right = scoreMatrix[(cur_row) * cols + cur_col - 1] + GAP;

	// calculate score from going diagonal
	int	diag = 0;
	if (seq1[cur_row-1] == seq2[cur_col-1]) {
		//seq1 row and seq2 col nucleotides match
		diag = scoreMatrix[(cur_row -1) * cols + (cur_col - 1)] + MATCH;
	}
	else {
		//seq1 row and seq2 col nucleotides do not match
		diag = scoreMatrix[(cur_row-1)* cols + cur_col-1] + MISMATCH;
	}

	// update cell value scoreand trace values based on highest score
	if (diag >= down && diag >= right) {
		//diag produced best score
		scoreMatrix[cur_row *cols + cur_col] = diag;
	}
	else if (right >= diag && right >= down) {
		// right produced best score
		scoreMatrix[cur_row *cols + cur_col] = right;
	}
	else {
		// down produced best score
		scoreMatrix[cur_row *cols + cur_col] = down;
	}
}

/*
	Launch a single block of 32 threads. 
	Iterates through first & last 32 diagonals of the matrix 
	Single warp, no __syncthreads needed.
*/
__global__ void warp_level_diag(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols, bool first_half) {
	int cur_row, cur_col;
	int diagIndex;
	int maxNumDiagonals = rows; // Assuming a square matrix . TODO: make this work for rectangular matrices too
	for (int i =0; i < 32; i++) {
		// For the first 32 diagonals of the matrix
		// Calculate the current row & col for each iteration
		if (first_half) {
			diagIndex = i;
			cur_row = threadIdx.x+1;
			cur_col =  diagIndex - threadIdx.x;
		} 
		// For the last 32 diagonals of the matrix
		// Calculate the current row & col for each iteration
		else {
			diagIndex = 32 - i; 
			cur_row  = maxNumDiagonals - diagIndex + threadIdx.x;
			cur_col = cols -1  - threadIdx.x;
		}
		// If thread maps to one of the diagonal cells, compute the value 
		if(threadIdx.x < diagIndex) {
			compute_cell(scoreMatrix, seq1, seq2, rows, cols, cur_row, cur_col);
		}
	}
}

/*
	Iterates through first & last 1024 diagonals of the matrix 
	If the matrix width/lenth is smaller than 1024 then iterate 
	through the entire matrix.
	Threads need to be synched at each iteration of the diagonal.
*/
__global__ void block_level_diag(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols, bool first_half) {
	int cur_row, cur_col;
	
	int index = blockIdx.x *  blockDim.x * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;

	/*
		For the first 1024 diagonals of the matrix
		Calculate the current row & col for each iteration
		If the matrix width/lenth is smaller than 1024
		then iterate through the entire upper diagonal of the matrix
			Assuming a square matrix. TODO: make this work for rectangular matrices too
	*/
	int maxNumDiagonals = min(1024, rows-1);  
	 if (first_half) {
		for (int i =32; i < maxNumDiagonals; i++) {
			cur_row = index + 1;
			cur_col = i - index;		
			// If thread maps to one of the diagonal cells, compute the value 
			if(index < i){
				compute_cell(scoreMatrix, seq1, seq2, rows, cols, cur_row, cur_col);
			}
			__syncthreads();
		}
	} 
	
	/*
		For the last 1024 diagonals of the matrix,
		Calculate the current row & col for each iteration
		If the matrix width/lenth is smaller than 1024
		then iterate through the entire lower diagonal of the matrix
	*/
	else {
		for (int i =maxNumDiagonals-1; i >= 31; i--) {
			cur_row = (rows-1) - i + index ;
			cur_col = (rows-1) - index ;
			// If thread maps to one of the diagonal cells, compute the value 
			if(index <= i ) {
				compute_cell(scoreMatrix, seq1, seq2, rows, cols, cur_row, cur_col);			
			}
			__syncthreads();
		}
	}
}

/*
	Runs through a single digonal with multiple blocks, assigning each thread to compute a single cell.
*/
__global__ void multi_block_level_diag(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols, int diagIndex, bool first_half) {
	int cur_row, cur_col;
	int index = blockIdx.x *  blockDim.x * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;

	// The rest of the upper diagonal of the matrix
	if (first_half) {
		cur_row = index + 1;
		cur_col = diagIndex - index;

		// If thread maps to one of the diagonal cells, compute the value 
		if(index < diagIndex){
			compute_cell(scoreMatrix, seq1, seq2, rows, cols, cur_row, cur_col);
		}
	} 
	// The rest of the lower diagonal of the matrix. 
	else {
		cur_row = rows - 1 - diagIndex + index ;
		cur_col = rows - 1 - index ;

		// If thread maps to one of the diagonal cells, compute the value 
		if(index <= diagIndex ) {
			compute_cell(scoreMatrix, seq1, seq2, rows, cols, cur_row, cur_col);			
		}	
	}
}

void launch_multi_kernel_diag(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols) {
	initialize_matrix << < ceil(max(cols, rows) / (float)256), 256 >> > (scoreMatrix, rows, cols);
	dim3 dimBlock(32, 32, 1);
	dim3 dimGrid(ceil(rows / (double)1024), 1, 1);

	// Single Kernel launch to go trhough first 32 diagonals of the matrix
	warp_level_diag << < 1, 32 >> > (scoreMatrix, seq1, seq2, rows, cols, true);
	
	
	// Single Kernel launch to go trhough first 1024 diagonals of the matrix
	block_level_diag << < dimGrid, dimBlock >> > (scoreMatrix, seq1, seq2, rows, cols, true);

	// Rest of the diagonals in the upper corner of the matrix
	// Launches multiples block of 32x32 threads. Each Kernel call computes a diagonal.
	for(int i =1024; i < rows; i++) {
		multi_block_level_diag << < dimGrid, dimBlock >> > (scoreMatrix, seq1, seq2, rows, cols, i, true);
	}

	// Rest of the diagonals in the lower corner of the matrix
	// Launches multiples block of 32x32 threads. Each Kernel call computes a diagonal.
	for(int i = rows-2; i >= 1024; i--) {
		multi_block_level_diag << < dimGrid, dimBlock >> > (scoreMatrix, seq1, seq2, rows, cols, i, false);

	}
	
	// Single Kernel launch to go trhough last 1024 diagonals of the matrix
	block_level_diag << < dimGrid, dimBlock >> > (scoreMatrix, seq1, seq2, rows, cols, false);
	
	// Single Kernel launch to go trhough last 32 diagonals of the matrix
	warp_level_diag << < 1, 32 >> > (scoreMatrix, seq1, seq2, rows, cols, false);
}


/*
======================================
	Multi Kernel Implementation / End
======================================
*/


void computeOnDevice(int* scoreMatrix, char* seq1, char* seq2, int rows, int cols, int n, int kernel_mode) {

	// setup code
	int* d_scoreMatrix;
	char* d_seq1, * d_seq2;

	cudaMalloc((void**)&d_scoreMatrix, sizeof(int) * (rows * cols));
	cudaMalloc((void**)&d_seq1, sizeof(char) * (rows - 1));
	cudaMalloc((void**)&d_seq2, sizeof(char) * (cols - 1));
	cudaMemcpy(d_seq1, seq1, sizeof(char) * (rows - 1), cudaMemcpyHostToDevice);
	cudaMemcpy(d_seq2, seq2, sizeof(char) * (cols - 1), cudaMemcpyHostToDevice);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	printf("\tTiming 'Needleman_Wunsch_GPU' started\n");

	
	if(kernel_mode == 1 || kernel_mode == 2){
		// Launch kernel call for GPU calculation
		cudaEventRecord(start);
		for (int i = 0; i < n; i++) {
			if(kernel_mode==1){
				launch_diag_block_kernel(d_scoreMatrix, d_seq1, d_seq2, rows, cols);
			}
			else{
				launch_multi_kernel_diag(d_scoreMatrix, d_seq1, d_seq2, rows, cols);
			}
		}

		cudaEventRecord(stop);
		cudaEventSynchronize(stop);
	}

	else if(kernel_mode == 3){
		cudaEventRecord(start);

		launch_single_block_kernel(d_scoreMatrix, d_seq1, d_seq2, rows, cols);

		cudaEventRecord(stop);
		cudaEventSynchronize(stop);
	}

	printf("\tTiming 'Needleman_Wunsch_GPU' ended\n");

	float device_ms = 0;
	cudaEventElapsedTime(&device_ms, start, stop);

	printf("\t%i iterations = %f\n", n, device_ms / 1000);

	// Teardown code
	cudaMemcpy(scoreMatrix, d_scoreMatrix, sizeof(int) * (rows * cols), cudaMemcpyDeviceToHost);
	cudaFree(d_scoreMatrix);
	cudaFree(d_seq1);
	cudaFree(d_seq2);
}

#endif // #ifndef _NW_SINGLE_KERNEL_H_
