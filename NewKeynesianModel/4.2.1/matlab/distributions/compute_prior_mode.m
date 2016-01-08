function m = compute_prior_mode(hyperparameters,shape)
% This function computes the mode of the prior distribution given the (two, three or four) hyperparameters
% of the prior distribution.
%    
% INPUTS 
%   hyperparameters     [double]    1*n vector of hyper parameters.
%   shape               [integer]   scalar specifying the prior shape:
%                                     shape=1 => Beta distribution,
%                                     shape=2 => Gamma distribution,
%                                     shape=3 => Gaussian distribution,
%                                     shape=4 => Inverse Gamma (type 1) distribution,
%                                     shape=5 => Uniform distribution,
%                                     shape=6 => Inverse Gamma (type 2) distribution.
%                                     
% OUTPUTS 
%   m       [double]    scalar or 2*1 vector, the prior mode.
%
% REMARKS 
% [1] The size of the vector of hyperparameters is 3 when the Gamma or Inverse Gamma is shifted and 4 when 
%     the support of the Beta distribution is not [0,1].      
% [2] The hyperparameters of the uniform distribution are the lower and upper bounds.    
% [3] The uniform distribution has an infinity of modes. In this case the function returns the prior mean.
% [4] For the beta distribution we can have 1, 2 or an infinity of modes.

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
m = NaN ;
switch shape
  case 1
    if (hyperparameters(1)>1 && hyperparameters(2)>1)
        m = (hyperparameters(1)-1)/(hyperparameters(1)+hyperparameters(2)-2) ;
    elseif (hyperparameters(1)<1 && hyperparameters(2)<1)
        m = [ 0 ; 1 ] ;
    elseif ( hyperparameters(1)<1 && hyperparameters(2)>1-eps ) || ( abs(hyperparameters(1)-1)<2*eps && hyperparameters(2)>1 )
        m = 0;
    elseif ( hyperparameters(1)>1 && hyperparameters(2)<1+eps ) || ( abs(hyperparameters(1)-1)<2*eps && hyperparameters(2)<1 )
        m = 1;
    elseif ( abs(hyperparameters(1)-1)<2*eps && abs(hyperparameters(2)-1)<2*eps )% Uniform distribution!
        m = .5 ;
    end
    if length(hyperparameters)==4
        m = m*(hyperparameters(4)-hyperparameters(3)) + hyperparameters(3) ;
    end
  case 2
    if hyperparameters(1)<1
        m = 0;
    else
        m = (hyperparameters(1)-1)*hyperparameters(2);
    end
    if length(hyperparameters)>2
        m = m + hyperparameters(3);
    end
  case 3
    m = hyperparameters(1);
  case 4
    % s  = hyperparameters(1)
    % nu = hyperparameters(2)
    m = 1/sqrt((hyperparameters(2)+1)/hyperparameters(1));%sqrt((hyperparameters(2)-1)/hyperparameters(1))
    if length(hyperparameters)>2
        m = m + hyperparameters(3);
    end
  case 5
    m = .5*(hyperparameters(2)-hyperparameters(1)) ;
  case 6
    % s  = hyperparameters(1)
    % nu = hyperparameters(2)
    m = hyperparameters(1)/(hyperparameters(2)+2) ;
    if length(hyperparameters)>2
        m = m + hyperparameters(3) ;
    end
  otherwise
    error('Unknown prior shape!')
end