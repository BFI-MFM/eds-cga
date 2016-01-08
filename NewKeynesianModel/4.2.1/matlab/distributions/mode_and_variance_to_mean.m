function [mu, parameters] = mode_and_variance_to_mean(m,s2,distribution,lower_bound,upper_bound)
% This function computes the mean of a distribution given the mode and variance of this distribution.
%
%  INPUTS 
%    m                [double]    scalar, mode of the distribution.
%    s2               [double]    scalar, variance of the distribution.
%    distribution     [integer]   scalar for the distribution shape 
%                                    1 gamma
%                                    2 inv-gamma-2
%                                    3 inv-gamma-1
%                                    4 beta    
%    lower_bound      [double]    scalar, lower bound of the random variable support (optional).
%    upper_bound      [double]    scalar, upper bound of the random variable support (optional).
%    
%  OUTPUT 
%    mu               [double]    scalar, mean of the distribution.
%    parameters       [double]    2*1 vector, parameters of the distribution.
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% Check input aruments. 
if ~(nargin==3 || nargin==5 || nargin==4 )
    error('mode_and_variance_to mean:: 3 or 5 input arguments are needed!')
end

% Set defaults bounds.
if nargin==3
    switch distribution
      case 1
        lower_bound = 0;
        upper_bound = Inf;
      case 3
        lower_bound = 0;
        upper_bound = Inf;
      case 2
        lower_bound = 0;
        upper_bound = Inf;
      case 4
        lower_bound = 0;
        upper_bound = 1;
      otherwise
        error('Unknown distribution!')
    end
end
if nargin==4
    switch distribution
      case 1
        upper_bound = Inf;
      case 3
        upper_bound = Inf;
      case 2
        upper_bound = Inf;
      case 4
        upper_bound = 1;
      otherwise
        error('Unknown distribution!')
    end
end


if (distribution==1)% Gamma distribution
    if m<lower_bound
        error('The mode has to be greater than the lower bound!')
    end
    if (m-lower_bound)<1e-12
        error('The gamma distribution should be specified with the mean and variance.')
    end        
    m = m - lower_bound ;
    beta  = -.5*m*(1-sqrt(1+4*s2/(m*m))) ;
    alpha = (m+beta)/beta ;
    parameters(1) = alpha;
    parameters(2) = beta;
    mu = alpha*beta + lower_bound ;
    return
end

if (distribution==2)% Inverse Gamma - 2 distribution
    if m<lower_bound+2*eps
        error('The mode has to be greater than the lower bound!')
    end
    m = m - lower_bound ;
    if isinf(s2)
        nu = 4;
        s  = 2*m;
    else
        delta = 2*(m*m/s2);
        poly = [ 1 , -(8+delta) , 20-4*delta , -(16+4*delta) ];
        all_roots  = roots(poly);
        real_roots = all_roots(find(abs(imag(all_roots))<2*eps));
        nu = real_roots(find(real_roots>2));
        s  = m*(nu+2);
    end
    parameters(1) = nu;
    parameters(2) = s;
    mu = s/(nu-2) + lower_bound;
    return
end

if (distribution==3)% Inverted Gamma 1 distribution
    if m<lower_bound+2*eps
        error('The mode has to be greater than the lower bound!')
    end
    m = m - lower_bound ;
    if isinf(s2)
        nu = 2;
        s  = 3*(m*m);
    else
        [mu, parameters] = mode_and_variance_to_mean(m,s2,2);
        nu = sqrt(parameters(1));
        nu2 = 2*nu;
        nu1 = 2;
        err = s2/(m*m) - (nu+1)/(nu-2) + .5*(nu+1)*(gamma((nu-1)/2)/gamma(nu/2))^2;
        while abs(nu2-nu1) > 1e-12
            if err < 0
                nu1 = nu;
                if nu < nu2
                    nu = nu2;
                else
                    nu = 2*nu;
                    nu2 = nu;
                end
            else
                nu2 = nu;
            end
            nu =  (nu1+nu2)/2;
            err = s2/(m*m) - (nu+1)/(nu-2) + .5*(nu+1)*(gamma((nu-1)/2)/gamma(nu/2))^2;
        end
        s = (nu+1)*(m*m) ;
    end
    parameters(1) = nu;
    parameters(2) = s;
    mu = sqrt(.5*s)*gamma(.5*(nu-1))/gamma(.5*nu) + lower_bound ;
    return
end

if (distribution==4)% Beta distribution
    if m<lower_bound
        error('The mode has to be greater than the lower bound!')
    end
    if m>upper_bound
        error('The mode has to be less than the upper bound!')
    end
    if (m-lower_bound)<1e-12
        error('The beta distribution should be specified with the mean and variance.')
    end
    if (upper_bound-m)<1e-12
        error('The beta distribution should be specified with the mean and variance.')
    end
    ll = upper_bound-lower_bound;
    m  = (m-lower_bound)/ll;
    s2 = s2/(ll*ll);
    delta = m^2/s2;
    poly = NaN(1,4);
    poly(1) = 1;
    poly(2) = 7*m-(1-m)*delta-3;
    poly(3) = 16*m^2-14*m+3-2*m*delta+delta;
    poly(4) = 12*m^3-16*m^2+7*m-1;
    all_roots = roots(poly);
    real_roots = all_roots(find(abs(imag(all_roots))<2*eps));
    idx = find(real_roots>1);
    if length(idx)>1
        error('Multiplicity of solutions for the beta distribution specification.')
    elseif isempty(idx)
        disp('No solution for the beta distribution specification. You should reduce the variance.')
        error();
    end
    alpha = real_roots(idx);
    beta = ((1-m)*alpha+2*m-1)/m;
    parameters(1) = alpha;
    parameters(2) = beta;
    mu = alpha/(alpha+beta)*(upper_bound-lower_bound)+lower_bound;
    return
end