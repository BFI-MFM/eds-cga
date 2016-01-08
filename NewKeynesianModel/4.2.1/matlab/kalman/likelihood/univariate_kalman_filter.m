function [LIK, llik] = univariate_kalman_filter(T,R,Q,H,P,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations)
% Computes the likelihood of a stationnary state space model (univariate approach).
%
% INPUTS
%    T                            [double]    mm*mm transition matrix of the state equation.
%    R                            [double]    mm*rr matrix, mapping structural innovations to state variables.
%    Q                            [double]    rr*rr covariance matrix of the structural innovations.
%    H                            [double]    pp*1 (zeros(pp,1) if no measurement errors) variances of the measurement errors. 
%    P                            [double]    mm*mm variance-covariance matrix with stationary variables
%    Y                            [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                        [integer]   scalar, likelihood evaluation starts at 'start'.
%    mf                           [integer]   pp*1 vector of indices.
%    kalman_tol                   [double]    scalar, tolerance parameter (rcond).
%    riccati_tol                  [double]    scalar, tolerance parameter (riccati iteration).
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    number_of_observations       [integer]   scalar.
%    no_more_missing_observations [integer]   scalar.
%
% OUTPUTS
%    LIK        [double]    scalar, MINUS loglikelihood
%    llik        [double]    vector, density of observations in each period.
%
% REFERENCES
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series
%   Analysis, vol. 24(1), pp. 85-98).
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

% Copyright (C) 2004-2010 Dynare Team
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

pp     = size(Y,1);                            % Number of observed variables 
mm     = size(T,1);                            % Number of variables in the state vector.
smpl   = size(Y,2);                            % Number of periods in the dataset.
a      = zeros(mm,1);                          % Initial condition of the state vector.
QQ     = R*Q*transpose(R);
t      = 0;
lik    = zeros(smpl,1); 
llik    = zeros(smpl,pp); 
notsteady   = 1;
l2pi = log(2*pi);

while notsteady && t<smpl
    t  = t+1;
    d_index = data_index{t};
    MF = mf(d_index);
    oldP = P;
    for i=1:length(MF)
        prediction_error = Y(d_index(i),t) - a(MF(i));
        Fi = P(MF(i),MF(i)) + H(d_index(i));
        if Fi > kalman_tol
            Ki     = P(:,MF(i))/Fi;
            a      = a + Ki*prediction_error;
            P      = P - (Fi*Ki)*transpose(Ki);
            llik(t,i) = log(Fi) + prediction_error*prediction_error/Fi ...
                     + l2pi;
            lik(t) = lik(t) + llik(t,i);
        end
    end
    a = T*a;
    P = T*P*transpose(T) + QQ;
    if t>no_more_missing_observations
        notsteady = max(max(abs(P-oldP)))>riccati_tol;
    end
end

% Steady state kalman filter.
while t < smpl
    PP = P;
    t = t+1;
    for i=1:pp
        prediction_error = Y(i,t) - a(mf(i));
        Fi   = PP(mf(i),mf(i)) + H(i);
        if Fi > kalman_tol
            Ki = PP(:,mf(i))/Fi;
            a  = a + Ki*prediction_error;
            PP  = PP - (Fi*Ki)*transpose(Ki);
            llik(t,i) = log(Fi) + prediction_error*prediction_error/Fi ...
                     + l2pi;
            lik(t) = lik(t) + llik(t,i);
        end
    end
    a = T*a;
end

lik = lik/2;
llik = llik/2;

LIK = sum(lik(start:end));