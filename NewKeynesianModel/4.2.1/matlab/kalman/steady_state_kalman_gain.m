function [K,iF,P] = steady_state_kalman_gain(T,R,Q,H,mf)
% Given the invariant state space representation of a model, this
% function computes the gain matrix and the covariance matrix of the
% state vector at the steady state of the kalman filter. 
% 
% INPUTS 
%   T   [double]    m*m transition matrix of the state vector.  
%   R   [double]    m*q matrix (q is the number of structural innovations).
%   Q   [double]    q*q covariance matrix of the structural innovations.
%   H   [double]    p*p covariance matrix of the measurement error.
%   mf  [integer]   p*1 vector, indices for the observed variables
%    
% OUTPUTS 
%   K   [double]    kalman gain matrix.
%   P   [double]    covariance matrix of the state vector.
%               
% SPECIAL REQUIREMENTS
%   Needs a solver for Riccati equations (dare.m)

% Copyright (C) 2004-2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

m = length(T);
p = length(mf);
Z = build_selection_matrix(mf,m,p);

if isempty(H)
    H = zeros(p,p);
end

QQ = R*Q*transpose(R);

P  = dare(T,transpose(Z),QQ,H);

iF = inv(Z*P*transpose(Z)+H);
K  = T*P*transpose(Z)*iF;