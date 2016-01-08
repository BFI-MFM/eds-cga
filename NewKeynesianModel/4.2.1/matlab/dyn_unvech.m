function Matrix = dyn_unvech(Vector)
% This function implements the unvech operator.
% 
% INPUTS 
%   Vector             [double]   a m*1 vector.
%    
% OUTPUTS 
%   Matrix             [double]   a n*n symetric matrix, where n solves n*(n+1)/2=m.

% Copyright (C) 2010 Dynare Team
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

m = length(Vector);
n = (sqrt(1+8*m)-1)/2;
b = 0;
Matrix = NaN(n,n);
for col = 1:n
    idx = 1:col;
    Matrix(1:col,col) = Vector(b+idx);
    Matrix(col,1:col) = transpose(Matrix(1:col,col));
    b = b+length(idx);
end