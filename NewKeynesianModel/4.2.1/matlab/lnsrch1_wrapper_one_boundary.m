function r = lnsrch1_wrapper_one_boundary(ya, y_index, fname, y, x, params, steady_state, it_)
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
% OUTPUTS
%   r                   [vector]        The residuals of the current block      
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
y(it_, :) = ya;
[r, y, g1, g2, g3]=feval(fname, y, x, params, steady_state, it_, 0);
