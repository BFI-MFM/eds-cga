function B = rand_matrix_normal(n, p, M, Omega_lower_chol, Sigma_lower_chol)

% function B = rand_matrix_normal(n, p, M, Omega_lower_chol, Sigma_lower_chol)
% Pseudo random matrices drawn from a matrix-normal distribution
% B ~ MN_n*p(M, Omega, Sigma) 
% Equivalent to vec(B) ~ N(vec(Mu), kron(Omega, Sigma))
%
% INPUTS
%    n:                 row
%    p:                 column
%    M:                 (n*p) matrix, mean
%    Omega_lower_chol:  (p*p), lower Cholesky decomposition of Omega,
%                       (Omega_lower_chol = chol(Omega, 'lower'))
%    Sigma_lower_chol:  (n*n), lower Cholesky decomposition of Sigma,
%                       (Sigma_lower_chol = chol(Sigma, 'lower'))
%    
% OUTPUTS
%    B:                 (n*p) matrix drawn from a Matrix-normal distribution
%        
% SPECIAL REQUIREMENTS
%    Same notations than: http://en.wikipedia.org/wiki/Matrix_normal_distribution

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

B1 = randn(n * p, 1);
B2 = kron(Omega_lower_chol, Sigma_lower_chol) * B1;
B3 = reshape(B2, n, p);
B = B3 + M;
