function [err, D] = A_times_B_kronecker_C(A,B,C,fake)
%function [err, D] = A_times_B_kronecker_C(A,B,C)
% Computes A * kron(B,C). 
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

% Chek number of inputs and outputs.
if nargin>4 || nargin<3 || nargout~=2
    error('A_times_B_kronecker_C takes 3 or 4 input arguments and provides exactly 2 output arguments.')
end


% Get & check dimensions. Initialization of the output matrix.
[mA,nA] = size(A);
[mB,nB] = size(B);
if nargin == 4
    [mC,nC] = size(C);
    if mB*mC ~= nA
        error('Input dimension error!')
    end
    D = zeros(mA,nB*nC);
    loop = (mB*nB*mC*nC > 1e7);
else
    if mB*mB ~= nA
        error('Input dimension error!')
    end
    D = zeros(mA,nB*nB);
    loop = (mB*nB*mB*nB > 1e7);
end
% Computational part.
if loop
    if nargin == 4
        k1 = 1; 
        for i1=1:nB
            for i2=1:nC
                D(:,k1) = A * kron(B(:,i1),C(:,i2));
                k1 = k1 + 1; 
            end
        end
    else
        k1 = 1;
        for i1=1:nB
            for i2=1:nB
                D(:,k1) = A * kron(B(:,i1),B(:,i2));
                k1 = k1 + 1; 
            end
        end
    end
else
    if nargin == 4
        D = A * kron(B,C);
    else
        D = A * kron(B,B);
    end
end
err = 0;