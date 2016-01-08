function optimal_bandwidth = mh_optimal_bandwidth(data,number_of_draws,bandwidth,kernel_function) 
% This function gives the optimal bandwidth parameter of a kernel estimator
% used to estimate a posterior univariate density from realisations of a 
% Metropolis-Hastings algorithm. 
%
% INPUTS:
%   data               [double]  Vector (number_of_draws*1) of draws.
%   number_of_draws    [integer] Scalar, number of draws.
%   bandwidth          [integer] Scalar equal to 0,-1 or -2.    
%                                bandwidth =  0 => Silverman [1986] rule of thumb.
%                                bandwidth = -1 => Sheather and Jones [1991].
%                                bandwidth = -2 => Local bandwith parameters.                              
%   kernel_function    [string]  Name of the kernel function: 'gaussian', 'triweight',
%                                'uniform', 'triangle', 'epanechnikov', 'quartic', 
%                                'triweight' and 'cosinus'
%
% OUTPUTS:
%   optimal_bandwidth: [double]  Scalar or vector, optimal window width.
%   
% SPECIAL REQUIREMENTS:
%   none.
%
% REFERENCES:
%   [1] M. Skold and G.O. Roberts [2003], "Density estimation for the Metropolis-Hastings algorithm". 
%   [2] Silverman [1986], "Density estimation for statistics and data analysis". 

% Copyright (C) 2004-2011 Dynare Team
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

%% Kernel specifications.
if strcmpi(kernel_function,'gaussian')
    % Kernel definition
    k    = inline('inv(sqrt(2*pi))*exp(-0.5*x.^2)');
    % Second derivate of the kernel function
    k2   = inline('inv(sqrt(2*pi))*(-exp(-0.5*x.^2)+(x.^2).*exp(-0.5*x.^2))');
    % Fourth derivate of the kernel function
    k4   = inline('inv(sqrt(2*pi))*(3*exp(-0.5*x.^2)-6*(x.^2).*exp(-0.5*x.^2)+(x.^4).*exp(-0.5*x.^2))');
    % Sixth derivate of the kernel function
    k6   = inline(['inv(sqrt(2*pi))*(-15*exp(-0.5*x.^2)+45*(x.^2).*exp(-' ...
                   '0.5*x.^2)-15*(x.^4).*exp(-0.5*x.^2)+(x.^6).*exp(-0.5*x.^2))']);
    mu02 = inv(2*sqrt(pi));
    mu21 = 1;
elseif strcmpi(kernel_function,'uniform') 
    k    = inline('0.5*(abs(x) <= 1)');
    mu02 = 0.5;
    mu21 = 1/3;
elseif strcmpi(kernel_function,'triangle') 
    k    = inline('(1-abs(x)).*(abs(x) <= 1)');
    mu02 = 2/3;
    mu21 = 1/6;
elseif strcmpi(kernel_function,'epanechnikov') 
    k    = inline('0.75*(1-x.^2).*(abs(x) <= 1)');
    mu02 = 3/5;
    mu21 = 1/5;
elseif strcmpi(kernel_function,'quartic') 
    k    = inline('0.9375*((1-x.^2).^2).*(abs(x) <= 1)');
    mu02 = 15/21;
    mu21 = 1/7;
elseif strcmpi(kernel_function,'triweight') 
    k    = inline('1.09375*((1-x.^2).^3).*(abs(x) <= 1)');
    k2   = inline('(105/4*(1-x.^2).*x.^2-105/16*(1-x.^2).^2).*(abs(x) <= 1)');
    k4   = inline('(-1575/4*x.^2+315/4).*(abs(x) <= 1)');
    k6   = inline('(-1575/2).*(abs(x) <= 1)');
    mu02 = 350/429;
    mu21 = 1/9;
elseif strcmpi(kernel_function,'cosinus') 
    k    = inline('(pi/4)*cos((pi/2)*x).*(abs(x) <= 1)');
    k2   = inline('(-1/16*cos(pi*x/2)*pi^3).*(abs(x) <= 1)');
    k4   = inline('(1/64*cos(pi*x/2)*pi^5).*(abs(x) <= 1)');
    k6   = inline('(-1/256*cos(pi*x/2)*pi^7).*(abs(x) <= 1)');
    mu02 = (pi^2)/16;
    mu21 = (pi^2-8)/pi^2;
else
    disp('mh_optimal_bandwidth:: ');
    error('This kernel function is not yet implemented in Dynare!');
end


%% Get the Skold and Roberts' correction.
if bandwidth==0 || bandwidth==-1
    correction = correction_for_repeated_draws(data,number_of_draws);
else
    correction = 0;
end
%% Compute the standard deviation of the draws.
sigma = std(data);
%% Optimal bandwidth parameter.
if bandwidth == 0;  % Rule of thumb bandwidth parameter (Silverman [1986].
    h = 2*sigma*(sqrt(pi)*mu02/(12*(mu21^2)*number_of_draws))^(1/5);
    h = h*correction^(1/5);
elseif bandwidth == -1; % Sheather and Jones [1991] plug-in estimation of the optimal bandwidth parameter. 
    if strcmp(kernel_function,'uniform')      || ... 
            strcmp(kernel_function,'triangle')     || ... 
            strcmp(kernel_function,'epanechnikov') || ... 
            strcmp(kernel_function,'quartic')
        error(['I can''t compute the optimal bandwidth with this kernel...' ...
               'Try the gaussian, triweight or cosinus kernels.']);
    end 
    Itilda4 = 8*7*6*5/(((2*sigma)^9)*sqrt(pi));
    g3      = abs(2*correction*k6(0)/(mu21*Itilda4*number_of_draws))^(1/9);
    Ihat3   = 0;
    for i=1:number_of_draws
        Ihat3 = Ihat3 + sum(k6((data(i,1)-data)/g3));
    end     
    Ihat3 = - Ihat3/((number_of_draws^2)*g3^7);
    g2    = abs(2*correction*k4(0)/(mu21*Ihat3*number_of_draws))^(1/7);
    Ihat2 = 0;
    for i=1:number_of_draws
        Ihat2 = Ihat2 + sum(k4((data(i)-data)/g2));
    end
    Ihat2 = Ihat2/((number_of_draws^2)*g2^5);
    h     = (correction*mu02/(number_of_draws*Ihat2*mu21^2))^(1/5); % equation (22) in Skold and Roberts [2003]. 
elseif bandwidth == -2;     % Bump killing... I compute local bandwith parameters in order to remove 
                            % spurious bumps introduced by long rejecting periods.   
    if strcmp(kernel_function,'uniform')      || ... 
            strcmp(kernel_function,'triangle')     || ... 
            strcmp(kernel_function,'epanechnikov') || ... 
            strcmp(kernel_function,'quartic')
        error(['I can''t compute the optimal bandwidth with this kernel...' ...
               'Try the gaussian, triweight or cosinus kernels.']);
    end
    T = zeros(n,1);
    for i=1:n
        j = i;
        while j<= n && (data(j,1)-data(i,1))<2*eps;
            j = j+1;
        end     
        T(i) = (j-i);
        correction = correction + 2*T(i) - 1;
    end
    correction = correction/number_of_draws;
    Itilda4 = 8*7*6*5/(((2*sigma)^9)*sqrt(pi));
    g3      = abs(2*correction*k6(0)/(mu21*Itilda4*correction))^(1/9);
    Ihat3   = 0;
    for i=1:number_of_draws
        Ihat3 = Ihat3 + sum(k6((data(i,1)-data)/g3));
    end
    Ihat3 = -Ihat3/((n^2)*g3^7);
    g2    = abs(2*correction*k4(0)/(mu21*Ihat3*n))^(1/7);
    Ihat2 = 0;
    for i=1:number_of_draws;
        Ihat2 = Ihat2 + sum(k4((data(i)-data)/g2));
    end     
    Ihat2 = Ihat2/((number_of_draws^2)*g2^5);
    h = ((2*T-1)*mu02/(number_of_draws*Ihat2*mu21^2)).^(1/5); % h is a column vector (local banwidth parameters). 
else
    disp('mh_optimal_bandwidth:: ');
    error('Parameter bandwidth must be equal to 0, -1 or -2.');
end

optimal_bandwidth = h;



function correction = correction_for_repeated_draws(draws,n)
correction = 0;
for i=1:n
    j = i;
    while j<=n && ( draws(j,1) - draws(i,1) )<2*eps; 
        j = j+1;
    end
    correction = correction + 2*(j-i) - 1;
end
correction = correction/n;