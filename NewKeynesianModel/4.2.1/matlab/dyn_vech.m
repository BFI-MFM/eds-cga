function Vector = dyn_vech(Matrix)
% This function implements the vech operator.
% 
% INPUTS 
%   Matrix             [double]   a squared n*n symetric matrix.
%    
% OUTPUTS 
%   Vector             [double]   a n*(n+1)/2 vector.

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

n = length(Matrix);
b = 0;
Vector = NaN(n*(n+1)/2,1);
for col = 1:n
    idx = transpose(1:col);
    Vector(b+idx) = Matrix((col-1)*n+idx);
    b = b+length(idx);
end