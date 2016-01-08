function [] = Tracing()
% DESCRIPTION
% This function is used to test the correct execution of a matlab section
% on remote machine.
% 
% If no error happen the function simply create a file.
%
% INPUTS
% ...
% 
% OUTPUTS
% ...
%        
% SPECIAL REQUIREMENTS
%    none

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

fid = fopen('Tracing.txt','w+');
fclose (fid);

exit
