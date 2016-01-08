function mexErrCheck(mexFunctionName, err)
% function mexErrCheck(mexFunctionName, err)
% this function halts processing if err is equal to 1.
%
% INPUTS
%   mexFunctionName [char]    Name of the mexFunction
%   err             [double]  error code returned from mexFunction
%
% OUTPUTS
%   none.
%
% ALGORITHM
%   ...
%
% SPECIAL REQUIREMENTS
%   none.
%

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

if ~ischar(mexFunctionName) || ~isscalar(err)
    error('The first argument must be a char and the second a scalar');
end

if err
    error(['Error encountered in: ' mexFunctionName '.']);
end