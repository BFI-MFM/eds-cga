function [xparams, logpost] = GetOneDraw(type)

% function [xparams, logpost] = GetOneDraw(type)
% draws one row from metropolis
%
% INPUTS
%    type:      posterior
%               prior
%        
% OUTPUTS
%    xparams:   vector of estimated parameters (drawn from posterior distribution)
%    logpost:   log of the posterior density relative to this row
%        
% SPECIAL REQUIREMENTS
%    none

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

switch type
  case 'posterior'
    [xparams, logpost] = metropolis_draw(0);
  case 'prior'
    xparams = prior_draw(0);
    logpost = evaluate_posterior_kernel(xparams');
end