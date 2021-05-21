# needleman-wunsch

#### Kia Vang, Jane Huynh, Jeanine Thao, Amir Jalili, Farris Al-Humayani

## Project Description
Global sequence alignment is a common technique in bioinformatics to determine an optimal alignment between two biological nucleotide sequences. The Needleman-Wunsch (NW) algorithm is a well known dynamic programming method for computing this optimal alignment. However, computing large input sequences is time-intensive for a CPU-bound implementation of NW. In this group project, we have developed a GPU-based approach that shows significant speedup for input sequences of ~45,000 nucleotides, but with potential to align much larger sizes on a high end GPU. After testing with different sizes, our observations proved that our implementation can achieve a speedup up to ~48x for an input sequence of 37,000 with a block size of 32x32.


## How to Run Code
  $ make
  
  $ ./needleman-wunsch <value> 

## Resources

https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm

https://ieeexplore.ieee.org/document/7346733
