function B = dynare_squeeze(A);
% Same as matlab's squeeze function except that it also affects 2D arrays.

% Copyright (C) 2009 Dynare Team
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

sizA = size(A); 
dimA = length(sizA);
switch dimA
  case 1
    B = A;
  case 2
    if sizA(1)==1
        B = transpose(A);
    elseif sizA(2)==1
        B = A(:,1);
    else
        B = A;
    end
  otherwise
    B =  squeeze(A);
end