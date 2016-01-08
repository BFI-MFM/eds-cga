function aa = isconst(y)
% Returns 1 if vector y is constant, 0 otherwise.
%  
% INPUTS:
%   yy        [double]    n*1 vector.
%
% OUTPUTS
%   aa        [integer]   scalar equal to 1 or 0.
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2008-2009 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>
aa = 0;
if all(abs(y(2:end)-y(1:end-1))<1e-10)
    aa = 1;
end