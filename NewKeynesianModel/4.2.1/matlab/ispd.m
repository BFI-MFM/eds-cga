function test = ispd(A)

% function test = ispd(A)
% Tests if a square matrix is positive definite. 
% 
% INPUTS 
%   o A       [double]   a square matrix. 
%  
% OUTPUTS 
%   o test    [integer]  is equal to one if A is pd, 0 otherwise. 
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2007-2009 Dynare Team
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

m = length(A);% I do not test for a square matrix...
test = 1;

for i=1:m
    if ( det( A(1:i, 1:i) ) < 2.0*eps )
        test = 0;
        break
    end
end