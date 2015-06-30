#Kohonen Self-Similarity Map
Python implementation of the Kohonen Self-Similarity map.

At the moment, things are pretty barebones; you can only do 2-d self-similarity maps, and everything is in Python (C subroutines are in the process of being written). Some of the goals are to:

1. Parallelize training and transformation with both OpenMP and MPI
2. Do the heavy lifting with C
3. Provide several options to customize speed vs. accuracy (i.e. heirarchical grid search to take N^2 -> NlogN)
4. Maybe throw in some quick visualization functions for fun.

