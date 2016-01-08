function [LIK, llik] = univariate_diffuse_kalman_filter(T,R,Q,H,Pinf,Pstar,Y,start,Z,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations)
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
%    Z                            [double]    pp*mm, selection matrix or pp independant linear combinations.
%    kalman_tol                   [double]    scalar, tolerance parameter (rcond).
%    riccati_tol                  [double]    scalar, tolerance parameter (riccati iteration).
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    number_of_observations       [integer]   scalar.
%    no_more_missing_observations [integer]   scalar.
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
global options_

[pp ,smpl] = size(Y); 
mm = size(T,1);
a  = zeros(mm,1);
QQ = R*Q*transpose(R);
t  = 0;
lik = zeros(smpl,1);
llik=zeros(smpl,pp);
notsteady = 1;
crit = 1.e-6;
newRank = rank(Pinf,crit);
l2pi = log(2*pi);

while newRank && (t<smpl)
    t = t+1;
    d_index = data_index{t};
    Za = Z(d_index,:);
    for i=1:length(d_index)
        Zi = Z(d_index(i),:);
        prediction_error = Y(d_index(i),t) - Zi*a;
        Fstar = Zi*Pstar*Zi' + H(i);
        Finf  = Zi*Pinf*Zi';
        Kstar = Pstar*Zi';
        if Finf>kalman_tol && newRank
            Kinf   = Pinf*Zi';
            Kinf_Finf = Kinf/Finf;
            a      = a + Kinf_Finf*prediction_error;
            Pstar  = Pstar + Kinf*(Kinf_Finf'*(Fstar/Finf)) - Kstar*Kinf_Finf' ...
                     - Kinf_Finf*Kstar';
            Pinf   = Pinf - Kinf*Kinf_Finf';
            llik(t,i) = log(Finf) + l2pi;
            lik(t) = lik(t) + llik(t,i);
        elseif Fstar>kalman_tol
            llik(t,i) = log(Fstar) + prediction_error* ...
                     prediction_error/Fstar + l2pi;
            lik(t) = lik(t) + llik(t,i);
            a = a + Kstar*(prediction_error/Fstar);
            Pstar = Pstar - Kstar*(Kstar'/Fstar);
        end
    end
    if newRank
        oldRank = rank(Pinf,crit);
    else
        oldRank = 0;
    end
    a     = T*a;
    Pstar = T*Pstar*T'+QQ;
    Pinf  = T*Pinf*T';
    if newRank
        newRank = rank(Pinf,crit);
    end
    if oldRank ~= newRank
        disp('univariate_diffuse_kalman_filter:: T does influence the rank of Pinf!')   
    end
end

if (t==smpl)
    error(['univariate_diffuse_kalman_filter:: There isn''t enough information to estimate the initial conditions of the nonstationary variables']);
    LIK = NaN;
    return
end

while notsteady && (t<smpl)
    t = t+1;
    oldP = Pstar;
    d_index = data_index{t};
    for i=1:length(d_index)
        Zi = Z(d_index(i),:);
        prediction_error = Y(d_index(i),t) - Zi*a;
        Ki   = Pstar*Zi';
        Fi   = Zi*Ki + H(i);
        if Fi > kalman_tol
            a      = a + Ki*(prediction_error/Fi);
            Pstar  = Pstar - Ki*(Ki'/Fi);
            llik(t,i) = log(Fi) + prediction_error*prediction_error/Fi ...
                     + l2pi;
            lik(t) = lik(t) + llik(t,i);
        end
    end 
    a     = T*a;
    Pstar = T*Pstar*T' + QQ;
    if t>no_more_missing_observations
        notsteady = max(max(abs(Pstar-oldP)))>riccati_tol;
    end
end

while t < smpl
    t = t+1;
    Pstar = oldP;
    for i=1:pp
        Zi = Z(i,:);
        prediction_error = Y(i,t) - Zi*a;
        Fi   = Zi*Pstar*Zi'+H(i);
        if Fi > kalman_tol
            Ki     = Pstar*Zi';
            a      = a + Ki*prediction_error/Fi;
            Pstar  = Pstar - Ki*Ki'/Fi;
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
