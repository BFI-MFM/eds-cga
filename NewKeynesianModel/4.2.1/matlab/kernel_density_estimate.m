function [abscissa,f] = kernel_density_estimate(data,number_of_grid_points,number_of_draws,bandwidth,kernel_function) 
% Estimates a continuous density. 
% 
% INPUTS
%   data                  [double]  Vector (number_of_draws*1) of draws.
%   number_of_grid_points [integer] Scalar, number of points where the density is estimated.
%                                   This (positive) integer must be a power of two.  
%   number_of_draws       [integer] Scalar, number of draws.
%   bandwidth             [double]  Real positive scalar.                                  
%   kernel_function       [string]  Name of the kernel function: 'gaussian', 'triweight',
%                                   'uniform', 'triangle', 'epanechnikov', 'quartic', 
%                                   'triweight' and 'cosinus'
%
% OUTPUTS
%    abscissa             [double] Vector (number_of_grid_points*1) of values on the abscissa axis.
%    f:                   [double] Vector (number_of_grid_points*1) of values on the ordinate axis, 
%                                  (density estimates).
%        
% SPECIAL REQUIREMENTS
%    none.
%
% REFERENCES
%    A kernel density estimator is used (see Silverman [1986], "Density estimation for statistics and data analysis")
%    The code is adapted from Anders Holtsberg's matlab toolbox (stixbox). 
%

% Copyright (C) 2004-2008 Dynare Team
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

if min(size(data))>1
    error('kernel_density_estimate:: data must be a one dimensional array.');
else
    data = data(:);
end

test = log(number_of_grid_points)/log(2);
if (abs(test-round(test)) > 1e-12)
    error('kernel_density_estimate:: The number of grid points must be a power of 2.');
end

%% Kernel specification.
if strcmpi(kernel_function,'gaussian') 
    kernel = @(x) inv(sqrt(2*pi))*exp(-0.5*x.^2);
elseif strcmpi(kernel_function,'uniform') 
    kernel = @(x) 0.5*(abs(x) <= 1); 
elseif strcmpi(kernel_function,'triangle') 
    kernel = @(x) (1-abs(x)).*(abs(x) <= 1);
elseif strcmpi(kernel_function,'epanechnikov') 
    kernel = @(x) 0.75*(1-x.^2).*(abs(x) <= 1);
elseif strcmpi(kernel_function,'quartic') 
    kernel = @(x) 0.9375*((1-x.^2).^2).*(abs(x) <= 1);
elseif strcmpi(kernel_function,'triweight') 
    kernel = @(x) 1.09375*((1-x.^2).^3).*(abs(x) <= 1);
elseif strcmpi(kernel_function,'cosinus') 
    kernel = @(x) (pi/4)*cos((pi/2)*x).*(abs(x) <= 1);
end

%% Non parametric estimation (Gaussian kernel should be used (FFT)).
lower_bound  = min(data) - (max(data)-min(data))/3;
upper_bound  = max(data) + (max(data)-min(data))/3;
abscissa = linspace(lower_bound,upper_bound,number_of_grid_points)';
inc = abscissa(2)-abscissa(1); 
xi  = zeros(number_of_grid_points,1);
xa  = (data-lower_bound)/(upper_bound-lower_bound)*number_of_grid_points; 
for i = 1:number_of_draws
    indx = floor(xa(i));
    temp = xa(i)-indx;
    xi(indx+[1 2]) = xi(indx+[1 2]) + [1-temp,temp]';
end
xk = [-number_of_grid_points:number_of_grid_points-1]'*inc;
kk = kernel(xk/bandwidth);
kk = kk / (sum(kk)*inc*number_of_draws);
f = ifft(fft(fftshift(kk)).*fft([ xi ; zeros(number_of_grid_points,1) ]));
f = real(f(1:number_of_grid_points));