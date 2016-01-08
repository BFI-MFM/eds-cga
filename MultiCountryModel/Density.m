% Density.m is a routine that estimates the density function in a given
% set of points by applying a multivariate kernel algorithm to the data;
% see formula (2) in Maliar and Maliar, (2015), "Merging simulation and 
% projection approaches to solve high-dimensional problems with an application
% to a new Keynesian model (henceforth, MM, 2015).
% -------------------------------------------------------------------------
% Inputs:  "Data"        is the matrix of d-dimensional data points; n-by-d
%          "Points"      is the matrix of data points in which the density   
%                        must be estimated; n_points-by-d; "Points" can be  
%                        a subset of "Data" or some other set of 
%                        d-dimensional points                    

% Output:  "Density"    is the estimated density function; n_points-by-1
%          "Points_PCn" is set of the normalized principal components of the 
%                        matrix of data for which the density must be estimated
%                       (needed for the EDS algorithm); n_points-by-d 
%          "Di_min "    is the distance from each point in "Data" to the 
%                       closest neighbors (needed for the EDS algorithm); 
%                       n_points-by-1
% -------------------------------------------------------------------------
% Copyright © 2015 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [Density, Points_PCn, Di_min] = Density(Data,Points) 
                               % "Data" is the matrix of all data; and "Points" 
                               % is the matrix of data points in which the   
                               % density is estimated ("Points" can be either  
                               % a subset of "Data" or some other set of  
                               % d-dimensional points)

[n,d] = size(Data);            % Infer the dimensionality of the data "Data"; 
                               % n-by-d

n_points = size(Points,1);     % Infer the number of rows in "Points"; the 
                               % number of columns is the same as in "Data", 
                               % i.e., d; n_points-by-d
   
Datan = (Data-ones(n,1)*mean(Data))./(ones(n,1)*std(Data)); 
                               % Transformation 1: normalize "Data" to zero 
                               % mean and unit standard deviation 

[U,S,V] = svd(Datan,0);        % Compute a singular-value decomposition (SVD)
                               % of matrix "Datan" using option "0", which 
                               % is "economy size" SVD in MATLAB; matrices 
                               % U, V and S are defined by Datan=U*S*V', 
                               % where U is n-by-d, S is d-by-d, V is d-by-d
                               
PC = Datan*V;                  % Transformation 2: compute principal components 
                               % using a linear change of variables
                               
PCn = PC./(ones(n,1)*std(PC)); % Transformation 3: normalize principal  
                               % components to unit standard deviation 
                               
Pointsn = (Points-ones(n_points,1)*mean(Data))./(ones(n_points,1)*std(Data)); 
                               % Transformation 1: normalize "CL" to zero 
                               % mean and unit standard deviation 
                               
Points_PC = Pointsn*V;         % Transformation 2: compute principal components 
                               % using a linear change of variables
                               
Points_PCn = Points_PC./(ones(n_points,1)*std(PC)); 
                               % Transformation 3: normalize principal  
                               % components to unit standard deviation 
                               
h_bar = n^(-1/(d+4));          % The bandwidth parameter; see formula (2) 
                               % in MM (2015)

Constant = 1/n/(2*pi)^(d/2)/h_bar^d;  
                               % A constant in formula (2) of MM (2015)
Density = zeros(n_points,1);   % Initialize the density; n_points-by-1
Di_min = zeros(n_points,1);    % Initialize the distance to the closest 
                               % point; n_points-by-1


for i = 1:n_points
    Di_2 = (ones(n,1)*Points_PCn(i,:)-PCn).^2*ones(d,1);
                               % Compute the squared distance of each point
                               % i from all n points; n-by-1
    Density(i,1) = Constant*ones(1,n)*exp(-0.5*Di_2/h_bar^2);
                               % Compute the density for the point i using
                               % formula (2) in MM (2015)
    % The distance to the closest neighbor is needed for the EDS algorithm 
    % --------------------------------------------------------------------
    Di_2(i,1) = inf;           % Set the squared distance between the point  
                               % i and itself to infinity to exclude it
                               % from the closest neighbor search 
    Di_min(i,1) = sqrt(min(Di_2));
                               % Find the distance from the point i to its
                               % closest neighbor 
end                                     