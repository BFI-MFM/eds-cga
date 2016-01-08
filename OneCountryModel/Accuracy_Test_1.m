% Accuracy_Test_1.m is a routine for evaluating accuracy of the one-agent 
% model's solutions, obtained by parameterizing capital policy function: it 
% computes approximation residuals in the optimality conditions on a given set 
% of points in the state space; see "Numerically Stable and Accurate 
% Stochastic Simulation Approaches for Solving Dynamic Economic Models" by 
% Kenneth L. Judd, Lilia Maliar and Serguei Maliar, (2011), Quantitative 
% Economics 2/2, 173–210 (henceforth, JMM, 2011).
%
% This version: July 14, 2011. First version: August 27, 2009.
% -------------------------------------------------------------------------
% Inputs:    "k" and "a" are, respectively, current-period capital and 
%            productivity levels, in the given set of points on which the 
%            accuracy is tested; 
%            "bk" are the coefficients of the capital policy function;
%            "IM" is the integration method for evaluating accuracy, 
%            IM=1,2,..,10=Gauss-Hermite quadrature rules with 1,2,...,10 
%            nodes in one dimension, respectively;
%            "PF" is the polynomial family, 0=Ordinary, 1=Hermite;
%            "zb" is a matrix of means and standard deviations of the state
%            variables, k and a; it is used to normalize these variables in 
%            the Hermite polynomial;            
%            "sigma", "rho", "beta", "gam", "alpha", "delta" and "A" are the 
%            parameters of the model;
%            "discard" is the number of data points to discard 
%
% Outputs:   "Residuals_mean" and "Residuals_max" are, respectively, the mean and
%            maximum absolute Euler equation residuals (in log10)
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [Residuals_mean Residuals_max time_test]  = Accuracy_Test_1(sigma,rho,beta,gam,alpha,delta,A,k,a,bk,D,IM,discard)

tic              % Start counting time needed to run the test

[n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(IM,1,sigma^2);
                 % n_nodes is the number of integration nodes, epsi_nodes  
                 % are integration nodes, and weight_nodes are integration 
                 % weights for Gauss-Hermite quadrature integration rule
                 % with IM nodes; for a unidimensional integral, n_nodes  
                 % coincides with IM


    for p = 1:size(a,1);     % For each point on which the accuracy is
                             % evaluated, ... 
        
        % Variables in point p 
        % ------------------------        
        k0 = k(p,1);         % Capital of period t
        a0 = a(p,1);         % Productivity level of period t 
        
        % Capital and consumption choices at t
        % ------------------------------------                      
        k1(1,1) =  Ord_Polynomial_N([k0(1,1) a0(1,1)],D)*bk;
        % Compute capital of period t+1 (chosen at t) using the corresponding 
        % capital policy function; 1-by-1
        c0(1,1) = A*k0(1,1)^alpha*a0(1,1) - k1(1,1) + (1-delta)*k0(1,1);
        % Consumption of period t
        
        % Capital and consumption choices at t+1
        %---------------------------------------
        a1 = a0.^rho.*exp(epsi_nodes); 
        % Productivity levels of period t+1; n_nodes-by-1 
        
        k1_dupl = ones(n_nodes,1)*k1; 
        % Duplicate k1 n_nodes times to create a matrix with n_nodes identical
        % rows; n_nodes-by-1 

        X1 = Ord_Polynomial_N([k1_dupl a1(:,1)],D);
        % Form a complete polynomial of degree D (at t+1) in the given point 
        
        k2(:,1) =  X1*bk;
        % Compute capital of period t+2 (chosen at t+1) using the fifth-
        % degree capital policy function; n_nodes-by-1 
        
        c1(:,1) =  A*k1(1,1)^alpha*a1(:,1) - k2(:,1) + (1-delta)*k1(1,1);
        % Consumption of period t+1; n_nodes-by-1

        % Approximation residuals in point p
        %--------------------------------  
        Residuals(p,1) = weight_nodes'*(beta*c1(1:n_nodes,1).^(-gam)./c0(1,1).^(-gam).*(1-delta+alpha*A*a1(1:n_nodes,1).*k1(1,1).^(alpha-1)))-1;
        % A unit-free Euler-equation approximation error
    end

 Residuals_mean = log10(mean(mean(abs(Residuals(1+discard:end,:)))));
 % Mean absolute Euler equation residuals (in log10)
 
 Residuals_max = log10(max(max(abs(Residuals(1+discard:end,:)))));    
 % Maximum absolute Euler equation residuals (in log10)
 

time_test = toc;   % Time needed to run the test  
