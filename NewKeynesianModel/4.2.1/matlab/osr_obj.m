function [loss,vx,info]=osr_obj(x,i_params,i_var,weights);
% objective function for optimal simple rules (OSR)

% Copyright (C) 2005-2009 Dynare Team
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

global M_ oo_ optimal_Q_ it_
%  global ys_ Sigma_e_ endo_nbr exo_nbr optimal_Q_ it_ ykmin_ options_

vx = [];
% set parameters of the policiy rule
M_.params(i_params) = x;

% don't change below until the part where the loss function is computed
it_ = M_.maximum_lag+1;
[dr,info] = resol(oo_.steady_state,0);

switch info(1)
  case 1
    loss = 1e8;
    return
  case 2
    loss = 1e8*min(1e3,info(2));
    return
  case 3
    loss = 1e8*min(1e3,info(2));
    return
  case 4
    loss = 1e8*min(1e3,info(2));
    return
  case 5
    loss = 1e8;
    return
  case 6
    loss = 1e8*min(1e3,info(2));
    return
  case 20
    loss = 1e8*min(1e3,info(2));
    return
  otherwise
end

vx = get_variance_of_endogenous_variables(dr,i_var);  
loss = weights(:)'*vx(:);















