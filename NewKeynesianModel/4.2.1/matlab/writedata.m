function writedata(fname)
% function writedata(fname)
% store endogenous and exogenous variables in a XLS spreadsheet file 
% INPUT
%   fname: name of the XLS file
% OUTPUT
%   none
% ALGORITHM
%   none
% SPECIAL REQUIREMENT
%   none

% Copyright (C) 2007-2009 Dynare Team
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

global M_ oo_

% xlswrite doesn't exist on Octave, and appeared in MATLAB 7.0
if exist('OCTAVE_VERSION') || matlab_ver_less_than('7.0')
    error('Function not supported on your version of MATLAB or Octave')
end

S=[fname '_endo.xls'];
n=size(oo_.endo_simul,2);
delete(S);
S=upper(cellstr(M_.endo_names));
S1=cellstr([num2str((1:n)')  char(65*ones(1,n))']);
xlswrite([fname '_endo'], S', 'endogenous', 'B1');
xlswrite([fname '_endo'], S1, 'endogenous', 'A2');
xlswrite([fname '_endo'], oo_.endo_simul', 'endogenous', 'B2');
S=[fname '_exo.xls'];
n=size(oo_.exo_simul,1);
delete(S);
S=upper(cellstr(M_.exo_names));
S1=cellstr([num2str((1:n)')  char(65*ones(1,n))']);
xlswrite([fname '_exo'], S','exogenous', 'B1');
xlswrite([fname '_exo'], S1, 'exogenous', 'A2');
xlswrite([fname '_exo'], oo_.exo_simul,'exogenous', 'B2');
