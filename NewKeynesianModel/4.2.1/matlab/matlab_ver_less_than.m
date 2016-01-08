function r = matlab_ver_less_than(verstr)
% function r = matlab_ver_less_than(verstr)
%
% Returns 1 if current Matlab version is strictly older than
% the one given in argument.
%
% It basically does the same job than verLessThan(), which is
% only available since Matlab 7.4.
%
% Note that this function will fail under Octave.
%
% INPUTS
%    verstr: a string of the format 'x.y' or 'x.y.z'
%    
% OUTPUTS
%    r: 0 or 1
%
% SPECIAL REQUIREMENTS
%    none

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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

ver_struct = ver('matlab');
cur_verstr = ver_struct.Version;

r = get_ver_numeric(cur_verstr) < get_ver_numeric(verstr);


function x = get_ver_numeric(verstr)
nums = sscanf(verstr, '%d.%d.%d')';
if length(nums) < 3
    nums(3) = 0;
end
x = nums * [1; 0.01; 0.0001 ];
