# Randomized-Time-Riemannian-Manifold-Hamiltonian-Monte-Carlo

This repository contains the C implementation of the simulations in "Randomized Time Riemannian Manifold Hamiltonian Monte Carlo" by Whalley, Paulin, and Leimkuhler.

The code here allows simulation of RT-RMHMC and RMHMC for all the examples presented in the paper including the data for the cosmological covariance estimation from Joachimi, B.: Non-linear shrinkage estimation of large-scale structure covariance. Monthly Notices of the Royal Astronomical Society: Letters 466(1), 83–87 (2017). 

The RMHMC_Bayesian.c and RT-RMHMC_Bayesian.c are set up for the IAC estimation for the upper quadrant, one needs to change the parameters at the top of the file to match those in the data extraction notebook to perform the Monte Carlo average. To put the data in the provided format for monte carlo averages one needs to unzip the provided data files of Joachimi, B.: Non-linear shrinkage estimation of large-scale structure covariance. Monthly Notices of the Royal Astronomical Society: Letters 466(1), 83–87 (2017) and run the provided jupyter notebook named Data Extraction with your choice of parameters. The jupyter notebook also provides a normalisation of the data.

Before compiling the code, one needs to insert their directory path into the RT-RMHMC_Bayesian.c and RMHMC_Bayesian.c files for importing the data. This is in place of "xxxxxx".

To run the C code you first need to compile by running
gcc Filenmae.c normal.c -lm -O3 -o out

To print the results to a file and run the code you can use
./out > dat.txt

If you are using the code please cite "Randomized Time Riemannian Manifold Hamiltonian Monte Carlo" by Whalley, Paulin, and Leimkuhler.
