function a=indnv(x,y)

% function a=indnv(x,y)
% Finds the elements indices of one vector in an other one
%
% INPUTS
%    x:         column vector
%    y:         column vector
%
% OUTPUTS
%    a:         vector of elements position of x in y
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2009 Dynare Team
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

a = zeros(size(x)) ;

for i = 1:size(x,1)
    j = find( x(i) == y );
    if isempty(j)
        a(i) = NaN;
    else
        a(i) = j;
    end
end



