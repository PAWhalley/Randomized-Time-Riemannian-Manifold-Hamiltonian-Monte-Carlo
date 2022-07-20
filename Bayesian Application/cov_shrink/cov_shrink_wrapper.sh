#!/bin/bash
# wrapper script for running NERCOME shrinkage covariance estimator
# Benjamin Joachimi, Nov 2016

### settings ###
nd=10           # dimension of input data vector
nr=100          # number of realisations of input data vector
nsplit=10       # number of split locations to test
nav=100         # number of samples to average over

ident=test      # file identifier
path=./         # path to files
param_file=./cov_shrink_${ident}.param  # parameter file name
data_file='-'  #data_${ident}.dat      # input data file; '-' for running with internally generated matrix with Gaussian iid variables
norm_file='-'   # file with normalisation; '-' for unit normalisation of covariance



### code ###
make cov_shrink
paramfile=${path}/${param_file}

# generate parameter file
echo $path > $paramfile
echo $data_file >> $paramfile
echo $ident >> $paramfile
echo $nd >> $paramfile
echo $nr >> $paramfile
echo $nsplit >> $paramfile
echo $nav >> $paramfile
echo $norm_file >> $paramfile


./cov_shrink $paramfile


# log settings
cat $0 | head -n 17 | tail -n 13 > ${path}/cov_shrink_${ident}.log
cat cov_shrink.c | head -n 26 | tail -n 7 >> ${path}/cov_shrink_${ident}.log
