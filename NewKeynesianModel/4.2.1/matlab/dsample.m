function dsample(s1,s2)
% function dsample(s1,s2)
% This optional command permits to reduce the number of periods considered in following output commands.
% If only one argument is provided, output is from period 1 to the period specified in the DSAMPLE command. 
% If two arguments are present output is done for the interval between the two periods.
% DSAMPLE without arguments reset the sample to the one specified by PERIODS
%
% INPUTS
%    s1:      first period
%    s2:      last period
%    
% OUTPUTS
%    none
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

global options_

options_.smpl = zeros(2,1) ;

if nargin == 0
    options_.smpl(1) = 1 ;
    options_.smpl(2) = options_.periods ;
elseif nargin == 1
    if s1 > options_.periods
        error('DSAMPLE: argument greater than number of periods');
    end
    options_.smpl(1) = 1 ;
    options_.smpl(2) = s1 ;
else
    if s1 > options_.periods || s2 > options_.periods
        error('DSAMPLE: one of the arguments is greater than number of periods');
    end
    options_.smpl(1) = s1 ;
    options_.smpl(2) = s2 ;
end

% 02/23/01 MJ added error checking