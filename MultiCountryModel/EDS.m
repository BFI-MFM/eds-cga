% EDS.m is a routine that constructs an epsilon distinguishable set 
% (EDS) for a given data set using Algorithm P_epsilon; see p. 10 of Maliar
% and Maliar (2015), "Merging simulation and projection approaches to solve
% high-dimensional problems with an application to a new Keynesian model 
% (henceforth, MM, 2015).
% -------------------------------------------------------------------------
% Inputs:  "Data" is the matrix of data on which the EDS set must be  
%                    constructed; n-by-d
%          "epsilon" is the minimum target distance between two points in 
%                    the EDS set
%          
% Output:  "EDS"     is the constructed EDS set on the original data; M-by-d
%          "EDS_PCn" is the constructed EDS set on principal components; M-by-d
% -------------------------------------------------------------------------
% Copyright © 2015 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% Nhe code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [EDS,EDS_PCn] = EDS(Data,epsilon) 
                               % "Data" is the matrix of data; "epsilon" is 
                               % the parameter that specifies the minimum 
                               % distance between any two data points in the
                               % EDS set
                               
[n,d] = size(Data);            % Infer the dimensionality of "Data"; n-by-d


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
                               
EDS_PCn = zeros(1,d);          % Initialize the epsilon distinguishable set
                               % of points; its largest possible size is 
                               % n-by-d
                               
eps2 = epsilon^2;              % Define the squared epsilon; we will compare 
                               % squared distances between points to the 
                               % squared epsilon (to save on cost of computing 
                               % the square root)

M = 0;                         % The number of points in EDS; initially, it 
                               % is zero
                               
N1 = n;                        % The number of data points left in PCn1 after  
                               % eliminating the points which are within the 
                               % the distance smaller than epsilon from the 
                               % EDS points; see Step 2 in Algorithm P_epsilon; 
                               % in the beginning, N1 is equal to the number 
                               % of data points, i.e., N1=n
                                                              
while N1>0;                    % While some data points are left unprocessed,... 
    
    M = M+1;                   % Increment the number of EDS points by 1
    EDS_PCn(M,:) = PCn(1,:);   % Add the first point of PCn to the EDS set  
    D_i2 = (ones(N1,1)*PCn(1,:)-PCn).^2*ones(d,1); 
                               % Compute the (squared) distance from the point 
                               % PCn(1,:) to all other points in PCn
    PCnsort = sortrows([D_i2, PCn],1); 
                               % Sort the points in PCn by the (squared) 
                               % distances to the given point PCn(1,:)
    nexc = ones(1,N1)*(D_i2<=eps2); 
                               % Select the points for which the (squared) 
                               % distance is inside (squared) epsilon and 
                               % compute the number of such points
    PCn = PCnsort(nexc+1:end,2:end);
                               % Eliminate the first "nexc" point from the 
                               % data since these points are within the 
                               % distance less than epsilon from the point 
                               % PCn(1,:)
    N1 = size(PCn,1);          % Compute the number of data points in the 
                               % data set after the elimination    
end

EDS_PC = EDS_PCn(1:M,:).*(ones(M,1)*std(PC)); 
                               % Backward transformation 3: re-normalize 
                               % the obtained EDS set "EDS_PCn" to have the  
                               % standard deviation std(PC) 

EDSn = EDS_PC*V';              % Backward transformation 2: use a linear 
                               % change of variables to express the EDS 
                               % set "EDS_PC" (obtained after the backward 
                               % transformation 3) in the same system of 
                               %  coordinates as the original normalized 
                               % variables

EDS = EDSn.*(ones(M,1)*std(Data))+(ones(M,1)*mean(Data)); 
                               % Backward transformation 1: re-normalize 
                               % the EDS set "EDSn" (obtained after the  
                               % backward transformation 2) to have the mean  
                               % and standard deviation as in the original 
                               % data, i.e., mean(Data) and std(Data),
                               % respectively