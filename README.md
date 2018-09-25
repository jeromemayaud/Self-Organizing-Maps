# Self-Organizing-Maps
A Matlab toolbox for Self-Organizing Maps (SOM) and others.

SOM Toolbox 2.0, a software library for Matlab 5 implementing the Self-Organizing Map algorithm is Copyright (C) 1999 by Esa Alhoniemi, Johan Himberg, Jukka Parviainen and Juha Vesanto.


# Running the SOM code
All files in the directory are needed to run the main Matlab file: 'data2kde2som'. This file converts categorical (i.e. binned) data into a kernel density estimation, and then runs this through SOM functionality by Vesanto et al. 

Two (CSV) inputs are needed to run 'data2kde2som' (1) bin_midpoints (the midpoints of each bin of your categorical data) and (2) the data that fits into your bins (each row representing data distributions for each point). 

# Performing Principal Component Anakysis (PCA)
There is also a file called 'pca_surrey' that performs PCA (to compare with SOM output) on the same two files.


