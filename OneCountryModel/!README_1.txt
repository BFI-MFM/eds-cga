The folder "OneCountryModel" contains MATLAB software that solves a 
one-country model using an epsilon distinguishable set algorithm (EDS) 
and cluster grid algorithm (CGA); it accompanies the article "Merging 
simulation and projection approaches to solve high-dimensional problems 
with an application to a new Keynesian model" by Lilia Maliar and Serguei 
Maliar (2015), Quantitative Economics 6/1, pages 1–47 (henceforth, MM, 2015). 

This software is based on that of Lilia  Maliar and Serguei Maliar for 
solving models using the generalized stochastic simulation algorithm (GSSA) 
method, as described in the paper "Numerically Stable and Accurate Stochastic
Simulation Approaches for Solving Dynamic Economic Models" by Kenneth L. 
Judd, Lilia Maliar and Serguei Maliar, (2011), Quantitative Economics 2/2, 
173–210 (henceforth, JMM, 2011). The modifications made are concerned with 
the construction of the grid.

This version: March 6, 2015. First version: April 17, 2010.
 
1. "Main_EDS_CGA_1.m"          is a main file for computing a solution to 
                               the one-country model using the EDS and 
                               CGA methods; 
2. "Accuracy_Test_1.m"         computes residuals of the optimality 
                               conditions of the one-country model on a   
                               given set of points in the state space; 
                               borrowed from JMM (2011);   
3. "Density.m"                 estimates the density function from a 
                               given set of points 
4. "Clusters.m"                constructs clusters from simulated series 
                               and computes clusters' centers (to be used 
                               as a grid) 
5. "EDS"                       constructs an epsilon distinguishable set 
                               for a given set of data (to be used as a 
                               grid)
6. "Ord_Polynomial_N.m"        constructs the sets of basis functions for                                 
                               ordinary polynomials of the degrees from 
                               one to five, for the N-country model; borrowed 
                               from JMM (2011)
7. "GH_Quadrature.m"           constructs integration nodes and weights for  
                               the Gauss-Hermite rules with the number of  
                               nodes in each dimension ranging from one to  
                               ten; borrowed from JMM (2011)                     
8. "epsi10200.mat"             contains the series of the productivity  
                               levels of length 10,200 that are used for 
                               computing solutions and for evaluating accuracy; 
                               borrowed from JMM (2011)

    To solve the model, execute "Main_EDS_CGA_1.m".

For updates and other related software, please, check the authors' web 
pages. For additional information, please, contact the corresponding 
author: Lilia Maliar, Department of Economics,  Stanford University, 
Stanford, CA 94305-6072, USA, maliarl@stanford.edu.

-------------------------------------------------------------------------
Copyright © 2015 by Lilia Maliar and Serguei Maliar. All rights reserved. 
The code may be used, modified and redistributed under the terms provided 
in the file "License_Agreement.txt".
-------------------------------------------------------------------------
