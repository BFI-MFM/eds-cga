function [LIK, llik] = ...
    univariate_kalman_filter_corr(T,R,Q,H,P,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations)
% Computes the likelihood of a stationnary state space model (univariate
% approach + correlated measurement errors).
%
% INPUTS
%    T                            [double]    mm*mm transition matrix of the state equation.
%    R                            [double]    mm*rr matrix, mapping structural innovations to state variables.
%    Q                            [double]    rr*rr covariance matrix of the structural innovations.
%    H                            [double]    pp*pp covariance matrix of the measurement error.  
%    P                            [double]    mm*mm variance-covariance matrix with stationary variables
%    Y                            [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                        [integer]   scalar, likelihood evaluation starts at 'start'.
%    Z                            [integer]   pp*mm selection matrix.
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

pp = size(Y,1);                            % Number of observed variables 
mm = size(T,1);                            % Number of variables in the state vector.
rr = size(R,2);                            % Number of structural innovations. 
smpl   = size(Y,2);                        % Number of periods in the dataset.
a      = zeros(mm+pp,1);                      % Initial condition of the state vector.
t      = 0;
lik    = zeros(smpl,1); 
llik    = zeros(smpl,pp); 
notsteady   = 1;

TT = zeros(mm+pp);
TT(1:mm,1:mm) = T;
QQ = zeros(rr+pp);
QQ(1:rr,1:rr) = Q;
QQ(rr+1:end,rr+1:end) = H;
RR = zeros(mm+pp,rr+pp);
RR(1:mm,1:rr) = R;
RR(mm+1:end,rr+1:end) = eye(pp);
PP  = zeros(mm+pp);
PP(1:mm,1:mm) = P;
PP(mm+1:end,mm+1:end) = H;
QQQQ = zeros(mm+pp);
RQR  = R*Q*R';
QQQQ(1:mm,1:mm) = RQR;
QQQQ(mm+1:end,mm+1:end) = H;
l2pi = log(2*pi);

while notsteady && t<smpl
    t  = t+1;
    d_index = data_index{t};
    MF = mf(d_index);
    oldPP = PP;
    for i=1:length(MF)
        prediction_error = Y(d_index(i),t) - a(MF(i)) - a( mm+i );
        Fi = PP(MF(i),MF(i)) + PP(mm+i,mm+i);         
        if Fi > kalman_tol
            llik(t,i) = log(Fi) + prediction_error*prediction_error/Fi ...
                     + l2pi;
            lik(t) = lik(t) + llik(t,i);
            Ki     = sum(PP(:,[MF(i) mm+i]),2)/Fi;
            a      = a + Ki*prediction_error;
            PP     = PP - (Ki*Fi)*transpose(Ki);
        end
    end
    a(1:mm) = T*a(1:mm);
    a(mm+1:end) = zeros(pp,1);
    PP(1:mm,1:mm) = T*PP(1:mm,1:mm)*transpose(T) + RQR;
    PP(mm+1:end,1:mm) = zeros(pp,mm);
    PP(1:mm,mm+1:end) = zeros(mm,pp);
    PP(mm+1:end,mm+1:end) = H;
    if t>no_more_missing_observations
        notsteady = max(max(abs(PP-oldPP)))>riccati_tol;
    end
end

% Steady state kalman filter.
while t < smpl
    PPPP = PP;
    t = t+1;
    for i=1:pp
        prediction_error = Y(i,t) - a(mf(i)) - a(mm+i);
        Fi   = PPPP(mf(i),mf(i)) + PPPP(mm+i,mm+i);
        if Fi > kalman_tol
            Ki = ( PPPP(:,mf(i)) + PPPP(:,mm+i) )/Fi;
            a  = a + Ki*prediction_error;
            PPPP  = PPPP - (Fi*Ki)*transpose(Ki);
            llik(t,i) = log(Fi) + prediction_error*prediction_error/Fi ...
                     + l2pi;
            lik(t) = lik(t) + llik(t,i);
        end
    end
    a(1:mm) = T*a(1:mm);
    a(mm+1:end) = zeros(pp,1);
end

lik = lik/2;
llik = llik/2;

LIK = sum(lik(start:end));