# Needleman-Wunsch 
#### EE5351 Applied Parallel Programming

#### Members: Kiana Vang, Jane Huynh, Jeanine Thao, Amir Jalili, Farris Al-Humayani


## Project Description
Global sequence alignment is a common technique in bioinformatics to determine an optimal alignment between two biological nucleotide sequences. The Needleman-Wunsch (NW) algorithm is a well known dynamic programming method for computing this optimal alignment. However, computing large input sequences is time-intensive for a CPU-bound implementation of NW. In this group project, we developed a GPU-based approach that shows significant speedup for input sequences of ~45,000 nucleotides, but with potential to align much larger sizes on a high end GPU. After testing with different sizes, our observations proved that our implementation can achieve a speedup up to ~48x for an input sequence of 37,000 with a block size of 32x32.

#### How to Run Code
Build the program:

    $ make

As a quick test, run the following command, which will compute two nucleotide sequences of length 5000 with the diagonal-block kernel:

    $ ./needleman_wunsch 1000 1

Options to run the program:

1. Run the program with 2 randomly generated nucleotide sequences with random lengths (1-100)

       $ ./needleman_wunsch
  
2. Run the program with 2 randomly generated nucleotide sequences that are specified sequence length

       $ ./needleman_wunsch [sequence size] [kernel mode]
  
3. Run the program with 2 nucleotide sequences contained in the specified text files

        $ ./needleman_wunsch [seq1.fna] [seq2.fna] [kernel mode]

The [kernel mode] parameter has 3 possible options, each of which will perform a different GPU implementation with varying speedups: 

1 = Multiple Kernels, Diagonal Blocks

2 = Multiple Kernels, Diagonal

3 = Single Kernel

## CPU Implementation

A dynamic programming implementation. Initially, the first column and row of the score matrix is filled by adding GAPs (-2) to each cell iteration. The algorithm then iterates through every row of the score matrix, with each thread calculating the score of an individual cell by accounting for the scores of the cells to the left, above, and to the top-left diagonal of the current cell.

## GPU Implementation

Multiple Kernel, Diagonal: Runs through a single digonal with multiple blocks, assigning each thread to compute a single cell.

Multiple Kernels, Diagonal Blocks: Runs through block-sized diagonals.

Single Kernel: Runs through the score matrix row by row in blocks. All cells in each block are computed through diagonals.

## Results

The following is an example of invoking the executable on the command line and the resulting status.

$ ./needleman_wunsch 5000 2<br />
    Timing 'Needleman_Wunsch_CPU' started<br />
        GetTimeOfDay Time (for 10 iterations) = 8.553<br />
        Clock Time        (for 10 iterations) = 8.54<br />
    Timing 'Needleman_Wunsch_CPU' ended<br />
        Timing 'Needleman_Wunsch_GPU' started<br />
        Timing 'Needleman_Wunsch_GPU' ended<br />
        10 iterations = 0.256918<br /><br />
    Test PASSED

## Resources

- [Needleman-Wunsch Wikipedia](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)
