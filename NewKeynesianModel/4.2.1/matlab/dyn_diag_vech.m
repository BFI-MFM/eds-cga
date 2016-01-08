function d = dyn_diag_vech(Vector)
% This function returns the diagonal elements of a symmetric matrix
% stored in vech form
% 
% INPUTS 
%   Vector             [double]   a m*1 vector.
%    
% OUTPUTS 
%   d                  [double]   a n*1 vector, where n solves n*(n+1)/2=m.

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
k = cumsum(1:n);
d = Vector(k);
