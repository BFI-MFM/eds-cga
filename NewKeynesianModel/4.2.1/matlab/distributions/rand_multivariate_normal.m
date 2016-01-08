function draw = rand_multivariate_normal(Mean,Sigma_upper_chol,n)
% Pseudo random draws from a multivariate normal distribution,
% \mathcal N_n(Mean,Sigma), with expectation Mean and variance Sigma.
%
% INPUTS 
%
%    Mean               [double]    1*n vector, expectation of the multivariate random variable.
%    Sigma_upper_chol   [double]    n*n matrix, upper triangular Cholesky decomposition of Sigma (the covariance matrix).
%    n                  [integer]   dimension.
%    
% OUTPUTS 
%    draw               [double]    1*n vector drawn from a multivariate normal distribution with expectation Mean and
%                                   covariance Sigma 
%        
% SPECIAL REQUIREMENTS

% Copyright (C) 2003-2009 Dynare Team
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

draw = Mean + randn(1,n) * Sigma_upper_chol;
