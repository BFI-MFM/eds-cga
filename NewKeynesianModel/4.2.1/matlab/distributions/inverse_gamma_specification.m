function [s,nu] = inverse_gamma_specification(mu,sigma,type)

% function [s,nu] = inverse_gamma_specification(mu,sigma,type)
% Specification of the inverse Gamma function parameters
% X ~ IG(s,nu)
%
% INPUTS
%    mu:      expectation
%    sigma:   standard deviation 
%    type=1:  inverse Gamma 1 
%    type=2:  inverse Gamma 2 

% OUTPUTS
%    s:       shape parameter
%    nu:      scale parameter 
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2010 Dynare Team
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

sigma2 = sigma^2;
mu2 = mu^2;

if type == 2;       % Inverse Gamma 2   
    nu   = 2*(2+mu2/sigma2);
    s    = 2*mu*(1+mu2/sigma2);
elseif type == 1;   % Inverse Gamma 1 
    if sigma2 < Inf;
        nu = sqrt(2*(2+mu2/sigma2));
        nu2 = 2*nu;
        nu1 = 2;
        err = 2*mu2*gamma(nu/2)^2-(sigma2+mu2)*(nu-2)*gamma((nu-1)/2)^2;
        while abs(nu2-nu1) > 1e-12
            if err > 0
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
            err = 2*mu2*gamma(nu/2)^2-(sigma2+mu2)*(nu-2)*gamma((nu-1)/2)^2;
        end
        s = (sigma2+mu2)*(nu-2);
    else;
        nu  = 2;
        s   = 2*mu2/pi;
    end;   
else;
    s  = -1;
    nu = -1;
end;

% 01/18/2004 MJ replaced fsolve with secant
%               suppressed chck
%               changed order of output parameters