function y=ff1_(x)

% function y=ff1_(x)
% splits the input argument x into endogenous and exogenous variables and calls the 'static' function
%
% INPUTS
%    x:          argument splitted between endogenous and exogenous
%        
% OUTPUTS
%    y:         'static' function residuals
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2008 Dynare Team
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

global it_ M_ oo_

n1 = size(x,1) - M_.exo_nbr;
oo_.exo_simul(it_+M_.maximum_lag-M_.maximum_lag,:) = x(n1+1:end)';
fh = str2func([M_.fname '_static']);
y=feval(fh,x(1:n1),oo_.exo_simul, M_.params);



