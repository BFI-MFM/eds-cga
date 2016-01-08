function save_params_and_steady_state(filename)
% function save_params_and_steady_state(filename)
%
% For all parameters, endogenous and exogenous variables, stores
% their value in a text file.
% * for parameters, the value is taken from the last parameter
%   initialization
% * for exogenous, the value is taken from the last initval block
% * for endogenous, the value is taken from the last steady state
%   computation (or, if no steady state has been computed, from the
%   last initval block)
% Note that no variable type is stored in the file, so that the values
% can be reloaded (with load_params_and_steady_state) in a setup where
% the variable types are different.
%  
% INPUTS
%   filename:   where to store the saved values
%  
% OUTPUTS
%   none
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ oo_

fid = fopen(filename, 'w');
if fid < 0
    error([ 'SAVE_PARAMS_AND_STEADY_STATE: Can''t open ' filename ]);
end

for i = 1:M_.param_nbr
    fprintf(fid, '%s %.16g\n', M_.param_names(i,:), M_.params(i));
end

for i = 1:M_.endo_nbr
    fprintf(fid, '%s %.16g\n', M_.endo_names(i,:), oo_.steady_state(i));
end

for i = 1:M_.exo_nbr
    fprintf(fid, '%s %.16g\n', M_.exo_names(i,:), oo_.exo_steady_state(i));
end

for i = 1:M_.exo_det_nbr
    fprintf(fid, '%s %.16g\n', M_.exo_det_names(i,:), oo_.exo_det_steady_state(i));
end

fclose(fid);
