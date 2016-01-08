function [LIK, lik] = kalman_filter(T,R,Q,H,P,Y,start,mf,kalman_tol,riccati_tol)
% Computes the likelihood of a stationnary state space model.
%
% INPUTS
%    T                      [double]    mm*mm transition matrix of the state equation.
%    R                      [double]    mm*rr matrix, mapping structural innovations to state variables.
%    Q                      [double]    rr*rr covariance matrix of the structural innovations.
%    H                      [double]    pp*pp (or 1*1 =0 if no measurement error) covariance matrix of the measurement errors. 
%    P                      [double]    mm*mm variance-covariance matrix with stationary variables
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                  [integer]   scalar, likelihood evaluation starts at 'start'.
%    mf                     [integer]   pp*1 vector of indices.
%    kalman_tol             [double]    scalar, tolerance parameter (rcond).
%    riccati_tol            [double]    scalar, tolerance parameter (riccati iteration).
%
% OUTPUTS
%    LIK        [double]    scalar, MINUS loglikelihood
%    lik        [double]    vector, density of observations in each period.
%
% REFERENCES
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series
%   Analysis, vol. 24(1), pp. 85-98).
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

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

smpl = size(Y,2);                               % Sample size.
mm   = size(T,2);                               % Number of state variables.
pp   = size(Y,1);                               % Maximum number of observed variables.
a    = zeros(mm,1);                             % State vector.
dF   = 1;                                       % det(F).
QQ   = R*Q*transpose(R);                        % Variance of R times the vector of structural innovations.
t    = 0;                                       % Initialization of the time index.
lik  = zeros(smpl,1);                           % Initialization of the vector gathering the densities.
LIK  = Inf;                                     % Default value of the log likelihood.
oldK = Inf;
notsteady   = 1;                                % Steady state flag.
F_singular  = 1;

while notsteady && t<smpl
    t  = t+1;
    v  = Y(:,t)-a(mf);
    F  = P(mf,mf) + H;
    if rcond(F) < kalman_tol
        if ~all(abs(F(:))<kalman_tol)
            return
        else
            a = T*a;
            P = T*P*transpose(T)+QQ;
        end
    else
        F_singular = 0;
        dF     = det(F);
        iF     = inv(F);
        lik(t) = log(dF)+transpose(v)*iF*v;
        K      = P(:,mf)*iF;
        a      = T*(a+K*v);
        P      = T*(P-K*P(mf,:))*transpose(T)+QQ;
        notsteady = max(abs(K(:)-oldK)) > riccati_tol;
        oldK = K(:);
    end
end

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

if t < smpl
    t0 = t+1;
    while t < smpl
        t = t+1;
        v = Y(:,t)-a(mf);
        a = T*(a+K*v);
        lik(t) = transpose(v)*iF*v;
    end
    lik(t0:smpl) = lik(t0:smpl) + log(dF);
end    

% adding log-likelihhod constants
lik = (lik + pp*log(2*pi))/2;

LIK = sum(lik(start:end)); % Minus the log-likelihood.