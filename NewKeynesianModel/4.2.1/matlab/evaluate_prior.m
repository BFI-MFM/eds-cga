function [ldens,parameters] = evaluate_prior(parameters)
% Evaluate the prior density at parameters.
%
% INPUTS
%    o parameters  a string ('posterior mode','posterior mean','posterior median','prior mode','prior mean') or a vector of values for 
%                  the (estimated) parameters of the model.
%    
%    
% OUTPUTS
%    o ldens       [double]  value of the logged prior density.
%    o parameters  [double]  vector of values for the estimated parameters.
%    
% SPECIAL REQUIREMENTS
%    None
%
% REMARKS
% [1] This function cannot evaluate the prior density of a dsge-var model...

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

global bayestopt_

if nargin==0
    parameters = 'posterior mode';
end

if ischar(parameters)
    switch parameters
      case 'posterior mode'
        parameters = get_posterior_parameters('mode');
      case 'posterior mean'
        parameters = get_posterior_parameters('mean');
      case 'posterior median'
        parameters = get_posterior_parameters('median');
      case 'prior mode'
        parameters = bayestopt_.p5(:);
      case 'prior mean'
        parameters = bayestopt_.p1;
      otherwise
        disp('eval_prior:: If the input argument is a string, then it has to be equal to:')
        disp('                ''posterior mode'', ')
        disp('                ''posterior mean'', ')
        disp('                ''posterior median'', ')
        disp('                ''prior mode'' or')
        disp('                ''prior mean''.')
        error
    end
end
clear('priordens');
ldens = priordens(parameters, bayestopt_.pshape, bayestopt_.p6, bayestopt_.p7, bayestopt_.p3, bayestopt_.p4);