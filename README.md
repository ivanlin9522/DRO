# DRO
This project helps implement Prof. Henry Lam's paper, "Recovering Best Statistical Guarantees via the Empirical Divergence-based Distributionally Robust Optimization".
The code is written in MATLAB.
The main function is multifunction.m, it handles different cases:
KaFlag is used to choose among chi-square with DOF of 1, k-1, and q_n.
Discrete is used to choose between discrete case and continuous case.
Oneside is used to choose between one-sided confidence interval and two-sided confidence interval.
SampleSize is n, which can be 20,30, and so on.

The prerequisite is that you have to download MATLAB and Mathematica, that's why we have the following dll file and lib file.
mathrun.h, ml64i3.dll, ml64i3m.lib, math.c and math.mexw64.
In order to run second_dev_cov_int function, please follow the instruction and know how to call Mathematica in MATLAB.
