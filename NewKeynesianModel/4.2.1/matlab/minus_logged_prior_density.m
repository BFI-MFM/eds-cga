function [f,fake] = minus_logged_prior_density(xparams,pshape,p6,p7,p3,p4)
% Evaluates minus the logged prior density.
% 
% INPUTS 
%   xparams    [double]   vector of parameters.
%   pshape     [integer]  vector specifying prior densities shapes.
%   p6         [double]   vector, first hyperparameter.
%   p7         [double]   vector, second hyperparameter.
%   p3         [double]   vector, prior's lower bound.
%   p4         [double]   vector, prior's upper bound. 
%
% OUTPUTS 
%   f          [double]  value of minus the logged prior density.

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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
fake = 1;
f = - priordens(xparams,pshape,p6,p7,p3,p4);