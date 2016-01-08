function ra = lnsrch1_wrapper_two_boundaries(ya, fname, y, y_index, x, ...
                                             params, steady_state, periods, y_kmin, y_size)
% wrapper for solve_one_boundary m-file when it is used with a dynamic
% model
%
% INPUTS
%   ya                  [vector]        The endogenous of the current block
%   y_index             [vector of int] The index of the endogenous variables of
%                                       the block
%   fname               [string]        name of the file containing the block
%                                       to simulate
%   y                   [matrix]        All the endogenous variables of the model
%   x                   [matrix]        All the exogenous variables of the model
%   params              [vector]        All the parameters of the model
%   periods             [int]           The number of periods
%   y_kmin              [int]           The maximum number of lag on en endogenous variables
%   y_size              [int]           The number of endogenous variables
%                                       in the current block
% OUTPUTS
%   ra                  [vector]        The residuals of the current block      
%  
% ALGORITHM
%   none.
%    
% SPECIAL REQUIREMENTS
%   none.
%  

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
% along with Dynare.  If not, see <http://www.gnu.org/licen

%reshape the input arguments of the dynamic function
y(y_kmin+1:y_kmin+periods, y_index) = reshape(ya',length(y_index),periods)';
[r, y, g1, g2, g3, b]=feval(fname, y, x, params, steady_state, periods, 0, y_kmin, y_size);
ra = reshape(r(:, y_kmin+1:periods+y_kmin),periods*y_size, 1);
