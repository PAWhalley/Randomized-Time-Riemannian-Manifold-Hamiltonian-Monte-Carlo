# Makefile for shrinkage covariance estimation code

CC = gcc
LB = -lgsl -lgslcblas -lm
OFLAG = -O3 -Wall
IFLAG =	
LFLAG =	


cov_shrink: cov_shrink.c
	$(CC) cov_shrink.c $(LB) $(BJLIB) $(OFLAG) $(IFLAG) $(LFLAG) -o cov_shrink

