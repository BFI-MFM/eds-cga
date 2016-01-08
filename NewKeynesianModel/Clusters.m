% Clusters.m is a routine that constructs clusters and computes the clusters' 
% centers for a given data set; see Judd, Maliar and Maliar, (2010), 
% "A Cluster-Grid Projection Method: Solving Problems with High Dimensionality", 
% NBER Working Paper 15965 (henceforth, JMM, 2010). 
% -------------------------------------------------------------------------
% Inputs:  "Data" is the matrix of data on which clusters must be constructed; 
%          T-by-M;
%          "n_clus" is is the number of clusters to be constructed;
%          1<=n_clus<=T

% Output:  "CL" are the centers of the constructed clusters on original
%          data; n_clus-by-M
%          "CL_PCn" are the centers of the constructed clusters on principal
%          components; n_clus-by-M
%          "clusters" is a vector that specifies to which cluster each data
%          point and PC is assigned; T-by-1;
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [CL,CL_PCn,clusters] = Clusters(Data,n_clus) 
                               % n_clus is the number of clusters to be 
                               % constructed
                               

[T,M] = size(Data);            % Infer the dimensionality of "Data"; T-by-M

Datan = (Data-ones(T,1)*mean(Data))./(ones(T,1)*std(Data)); 
                               % Transformation 1: normalize "Data" to zero 
                               % mean and unit standard deviation 

[U,S,V] = svd(Datan,0);        % Compute a singular-value decomposition (SVD)
                               % of matrix "Datan" using option "0", which 
                               % is "economy size" SVD in MATLAB; matrices 
                               % U, V and S are defined by Datan=U*S*V', 
                               % where U is T-by-M, S is M-by-M, V is M-by-M
                               
PC = Datan*V;                  % Transformation 2: compute principal components 
                               % using a linear change of variables; see JMM (2009)
                               
PCn = PC./(ones(T,1)*std(PC)); % Transformation 3: normalize principal  
                               % components to unit standard deviation 

clusters = clusterdata(PCn,'distance','euclidean','maxclust',n_clus,'linkage','ward'); 
                               % Distinguish n_clus clusters using hierarchical
                               % clustering algorithm with Euclidean distance 
                               % and Ward linkage; see JMM (2010)

for i = 1:n_clus;
    cluster_i =  (clusters==i);% Create a dummy variable to be used for 
                               % selecting the data points in cluster i
    CL_PCn(i,1:M) = sum(PCn.*(cluster_i*ones(1,M)))'/sum(cluster_i); 
                               % Compute the center of cluster i
end

CL_PC = CL_PCn.*(ones(n_clus,1)*std(PC)); 
                               % Backward transformation 3: re-normalize 
                               % the clusters' centers CL_PCn to have the 
                               % standard deviation std(PC) 

CLn = CL_PC*V';                % Backward transformation 2: use a linear 
                               % change of variables to express the clusters' 
                               % centers CL_PC in the same system of coordinates 
                               % as the original normalized variables

CL = CLn.*(ones(n_clus,1)*std(Data))+(ones(n_clus,1)*mean(Data)); 
                               % Backward transformation 1: re-normalize 
                               % the clusters' centers CLn to have the mean and 
                               % standard deviation as in the original data,
                               % i.e., mean(Data) and std(Data), respectively