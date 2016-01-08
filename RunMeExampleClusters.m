% RunMeExampleClusters.m is a simple example illustrating the construction of 
% clusters using a hierarchical clustering algorithm with Euclidean distance  
% and Ward linkage, as described in the article "Merging simulation and  
% projection approaches to solve high-dimensional problems with an application 
% to a new Keynesian model" by Lilia Maliar and Serguei Maliar (2015), 
% Quantitative Economics 6/1, pages 1–47 (henceforth, MM, 2015). 
%
% This version: March 19, 2015. First version: June 27, 2011.
% 
% -------------------------------------------------------------------------
% Copyright © 2015 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

clc;
clear all;

% 1. Parameters of the model
% ---------------------------
alpha   = 0.36;                 % Capital share 
beta    = 0.99;                 % Discount factor
d       = 1;                    % Depreciation rate
sigma   = 0.01;                 % Standard deviation for log shock
rho     = 0.95;                 % Persistence of log technology shock
A       = (1/beta-(1-d))/alpha; % Technology level            

% 2. Simulate time-series solution to the neoclassical growth model
% -----------------------------------------------------------------
% Under the assumption of full depreciation of capital and logarithmic utility 
% function, this model admits a closed form solution

% 2.1 Allocate memory for the simulated series
% --------------------------------------------
T    = 1000000;             % Length of simulation
a    = ones(T,1);           % Productivity levels
k    = ones(T,1);           % Capital

% 2.2 Simulate capital and productivity series
% --------------------------------------------
epsi = randn(T,1)*sigma;    % Draw productivity shocks 
for t = 1:T-1; 
    k(t+1,1) = alpha*beta*A*a(t)*k(t)^alpha;      % Capital
    a(t+1,:) = a(t,:).^rho.*exp(epsi(t,:));     % Productivity
end;

% 3. Reduce the number of simulated points by selecting each 100th point 
% -----------------------------------------------------------------------
% This is a cheap method that allows to obtain simulated points that are 
% uncorrelated (independent) across time

N = 10000;         % Select N = 10,000 points out of T = 1,000,0000 points 
a1 = ones(N,1);    % Allocate memory for productivity
k1 = ones(N,1);    % Allocate memory for capital
for j = 1:10000
    a1(j,:) = a(j*100,1);   
                   % Select each 100th point
    k1(j,:) = k(j*100,1);   
end

Data = [k1 a1];    % Data that will be used to construct the EDS grid
 
% 4. Estimating density function on principal components (PCs)
% ------------------------------------------------------------
[density PCn]  = Density(Data,Data);    % Estimate the density function in 
                                        % each data point
Sorted_PCn=sortrows([PCn density],3);   % Sort principal components by
                                        % density function
Sorted_Data=sortrows([Data density],3); % Sort Data by density function  

cutoff = 500;                           % The number of lowest density 
                                        % points to be dropped
Data1 = Sorted_Data(1+cutoff:end,1:2);  % Data after removing the lowest 
                                        % density points                                                                              
PCn1 = Sorted_PCn(1+cutoff:end,1:2);    % Principal components after 
                                        % removing the lowest density
                                        % points

% 5. Construct clusters on principal components (PCs)
% ---------------------------------------------------
n_clus = 10;                                % Number of clusters to construct

[CL,CL_PCn,clusters] = Clusters(Data1,n_clus); 
                                            % Construct clusters; "CL" and 
                                            % "CL_PCn" contain clusters' 
                                            % centersin the original 
                                            % coordinates centers and in 
                                            % PCs and "clusters" contains 
                                            % number of cluster to which 
                                            % each point of "PCn1" is assigned
                                            
PCn1_Sorted = sortrows([PCn1 clusters],3);  % Sort PCn1 by clusters   

for i=1:N
    index(i) = sum(clusters==i);            % Compute the number of elements 
end                                         % in each cluster

cindex = cumsum(index);                     % Compute a cumulative number of
                                            % elements; this is needed to
                                            % draw the clusters in the
                                            % figure
% 6. Figures 
% ----------
subplot(3,2,1);
plot(k1(:,1),a1(:,1),'.','MarkerSize',2), xlabel('k_t'), ylabel('a_t'),axis([0.82 1.18 0.85 1.15]);
title('Figure 1a. Simulated points');

subplot(3,2,2);
plot(Sorted_PCn(:,1),Sorted_PCn(:,2),'.','MarkerSize',2), xlabel('PC^1_t'), ylabel('PC^2_t'),axis([-5.5 5.5 -5.5 5.5]);
title('Figure 1b. Principal components (PCs)');

subplot(3,2,3);
axis([-4 4 -4 4]);
plot(Sorted_PCn(cutoff+1:10000,1),Sorted_PCn(cutoff+1:10000,2),'.',Sorted_PCn(1:cutoff,1),Sorted_PCn(1:cutoff,2),'x','MarkerSize',4), xlabel('PC^1_t'), ylabel('PC^2_t'),axis([-5.5 5.5 -5.5 5.5]);
title('Figure 1c. Dropped off low density points');

subplot(3,2,4);
plot(PCn1_Sorted(1:cindex(1),1),PCn1_Sorted(1:cindex(1),2),'.',PCn1_Sorted(cindex(1)+1:cindex(2),1),PCn1_Sorted(cindex(1)+1:cindex(2),2),'.',PCn1_Sorted(cindex(2)+1:cindex(3),1),PCn1_Sorted(cindex(2)+1:cindex(3),2),'.',PCn1_Sorted(cindex(3)+1:cindex(4),1),PCn1_Sorted(cindex(3)+1:cindex(4),2),'.',PCn1_Sorted(cindex(4)+1:cindex(5),1),PCn1_Sorted(cindex(4)+1:cindex(5),2),'.',PCn1_Sorted(cindex(5)+1:cindex(6),1),PCn1_Sorted(cindex(5)+1:cindex(6),2),'.',PCn1_Sorted(cindex(6)+1:cindex(7),1),PCn1_Sorted(cindex(6)+1:cindex(7),2),'.',PCn1_Sorted(cindex(7)+1:cindex(8),1),PCn1_Sorted(cindex(7)+1:cindex(8),2),'.',PCn1_Sorted(cindex(8)+1:cindex(9),1),PCn1_Sorted(cindex(8)+1:cindex(9),2),'.',PCn1_Sorted(cindex(9)+1:cindex(10),1),PCn1_Sorted(cindex(9)+1:cindex(10),2),'.','MarkerSize',10), xlabel('PC^1_t'), ylabel('PC^2_t'),axis([-5.5 5.5 -5.5 5.5]);
title('Figure 1d. Clusters on PCs');

subplot(3,2,5);
plot(CL_PCn(:,1),CL_PCn(:,2),'o','MarkerSize',4'), xlabel('PC^1_t'), ylabel('PC^2_t'),axis([-5.5 5.5 -5.5 5.5]);
title('Figure 1e. Cluster centers on PCs');

subplot(3,2,6);
axis([-4 4 -4 4]);
plot(CL(:,1),CL(:,2),'o','MarkerSize',4), xlabel('k_t'), ylabel('a_t'),axis([0.82 1.18 0.85 1.15]);
title('Figure 1f. Cluster centers on original data');