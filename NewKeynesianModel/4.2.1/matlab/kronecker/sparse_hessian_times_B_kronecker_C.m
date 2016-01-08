function [err, D] = sparse_hessian_times_B_kronecker_C(varargin)
%function [err, D] = sparse_hessian_times_B_kronecker_C(A,B,C, fake)
% Computes A * kron(B,C) where A is a sparse matrix.
%
% INPUTS
%   A   [double] mA*nA matrix.
%   B   [double] mB*nB matrix.
%   C   [double] mC*nC matrix.
%  
% OUTPUTS
%   err [double] scalar: 1 indicates failure, 0 indicates success
%   D   [double] mA*(nC*nB) or mA*(nB*nB) matrix.
%  
% ALGORITHM
%   none.    
%
% SPECIAL REQUIREMENTS
%   none.

% Copyright (C) 1996-2010 Dynare Team
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
A = varargin{1};
B = varargin{2};
C = varargin{3};
fake = varargin{nargin};
if nargout~=2
    error('sparse_hessian_times_B_kronecker_C provides exactly 2 output arguments.')
end

switch nargin
  case 4
    [fake,D] = A_times_B_kronecker_C(A,B,C,fake);
  case 3
    [fake,D] = A_times_B_kronecker_C(A,B,C);
  otherwise
    error('Two or Three input arguments required!')
end
err = 0;