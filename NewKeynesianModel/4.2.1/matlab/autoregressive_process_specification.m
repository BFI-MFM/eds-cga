function [InnovationVariance,AutoregressiveParameters] = autoregressive_process_specification(Variance,Rho,p)
% This function computes the parameters of an AR(p) process from the variance and the autocorrelation function
% (the first p terms) of this process.
%
% INPUTS 
%  [1] Variance                 [double]  scalar, variance of the variable.
%  [2] Rho                      [double]  p*1 vector, the autocorelation function: \rho(1), \rho(2), ..., \rho(p).
%  [3] p                        [double]  scalar, the number of lags in the AR process.
%
% OUTPUTS 
%  [1] InnovationVariance       [double]  scalar, the variance of the innovation.
%  [2] AutoregressiveParameters [double]  p*1 vector of autoregressive parameters.
%
% NOTES 
%
% The AR(p) model for {y_t} is:
%   
%           y_t = \phi_1 * y_{t-1} +  \phi_2 * y_{t-2} + ... +  \phi_p * y_{t-p} + e_t    
%
% Let \gamma(0) and \rho(1), ..., \rho(2) be the variance and the autocorrelation function of {y_t}. This function
% compute the variance of {e_t} and the \phi_i (i=1,...,p) from the variance and the autocorrelation function of {y_t}. 
% We know that:
%    
%    \gamma(0) = \phi_1 \gamma(1) + ... + \phi_p \gamma(p) + \sigma^2
%
% where \sigma^2 is the variance of {e_t}. Equivalently we have:
%
%    \sigma^2 = \gamma(0) (1-\rho(1)\phi_1 - ... - \rho(p)\phi_p)     
%
% We also have for any integer  h>0:
% 
%    \rho(h) = \phi_1 \rho(h-1) + ... + \phi_p \rho(h-p)
%
% We can write the equations for \rho(1), ..., \rho(p) using matrices. Let R be the p*p autocorelation
% matrix and v be the p*1 vector gathering the first p terms of the autocorrelation function. We have: 
%
%           v = R*PHI
%    
% where PHI is a p*1 vector with the autoregressive parameters of the AR(p) process. We can recover the autoregressive
% parameters by inverting the autocorrelation matrix: PHI = inv(R)*v.
% 
% This function first computes the vector PHI by inverting R and computes the variance of the innovation by evaluating
%
%           \sigma^2 = \gamma(0)*(1-PHI'*v)

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
AutoregressiveParameters = NaN(p,1);
InnovationVariance = NaN;
switch p
  case 1
    AutoregressiveParameters = Rho(1);
  case 2
    tmp = (Rho(2)-1)/(Rho(1)*Rho(1)-1);
    AutoregressiveParameters(1) = Rho(1)*tmp;
    AutoregressiveParameters(2) = 1-tmp;
  case 3
    t1 = 1/(Rho(2)-2*Rho(1)*Rho(1)+1);
    t2 = (1.5*Rho(1)-2*Rho(1)*Rho(1)*Rho(1)+.5*Rho(3))*t1;
    t3 = .5*(Rho(1)- Rho(3))/(Rho(2)-1);
    AutoregressiveParameters(1) = t2-t3-Rho(1);
    AutoregressiveParameters(2) = (Rho(2)*Rho(2)-Rho(3)*Rho(1)-Rho(1)*Rho(1)+Rho(2))*t1 ;
    AutoregressiveParameters(3) = t3-Rho(1)+t2;
  otherwise
    AutocorrelationMatrix = eye(p);
    for i=1:p-1
        AutocorrelationMatrix = AutocorrelationMatrix + diag(Rho(i)*ones(p-i,1),i);
        AutocorrelationMatrix = AutocorrelationMatrix + diag(Rho(i)*ones(p-i,1),-i);
    end
    AutoregressiveParameters = AutocorrelationMatrix\Rho;
end
InnovationVariance = Variance * (1-AutoregressiveParameters'*Rho);