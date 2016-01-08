function lpkern = evaluate_posterior_kernel(parameters,llik)
% Evaluate the prior density at parameters.
%
% INPUTS
%    o parameters  a string ('posterior mode','posterior mean','posterior median','prior mode','prior mean') or a vector of values for 
%                  the (estimated) parameters of the model.
%    
%    
% OUTPUTS
%    o lpkern      [double]  value of the logged posterior kernel.
%    
% SPECIAL REQUIREMENTS
%    None
%
% REMARKS
% [1] This function cannot evaluate the prior density of a dsge-var model...
% [2] This function use persistent variables for the dataset and the description of the missing observations. Consequently, if this function 
%     is called more than once (by changing the value of parameters) the sample *must not* change.    

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

[ldens,parameters] = evaluate_prior(parameters);
if nargin==1
    llik = evaluate_likelihood(parameters);
end
lpkern = ldens+llik;