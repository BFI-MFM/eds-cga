function [LIK, lik] = missing_observations_diffuse_kalman_filter(T,R,Q,H,Pinf,Pstar,Y,start,Z,kalman_tol,riccati_tol,...
                                                  data_index,number_of_observations,no_more_missing_observations)
% Computes the diffuse likelihood of a state space model.
%
% INPUTS
%    T           [double]      mm*mm transition matrix
%    R           [double]      mm*rr matrix
%    Q           [double]      rr*rr covariance matrix of the structural innovations.
%    H           [double]      pp*pp covariance matrix of the measurement errors (if H is equal to zero (scalar) there is no measurement error). 
%    Pinf        [double]      mm*mm matrix used to initialize the covariance matrix of the state vector.
%    Pstar       [double]      mm*mm matrix used to initialize the covariance matrix of the state vector.
%    Y           [double]      pp*smpl matrix of (detrended) data, where pp is the number of observed variables.
%    start       [integer]     scalar, likelihood evaluation starts at 'start'.
%    Z           [double]      pp*mm matrix, selection matrix or pp linear independant combinations of the state vector.
%    kalman_tol  [double]      scalar, tolerance parameter (rcond).
%    riccati_tol [double]      scalar, tolerance parameter (riccati iteration).
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    number_of_observations       [integer]   scalar.
%    no_more_missing_observations [integer]   scalar.
%             
% OUTPUTS
%    LIK:    MINUS loglikelihood
%    lik:    density vector in each period
%        
% REFERENCES
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series 
%   Analysis, vol. 24(1), pp. 85-98). 

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

[pp,smpl] = size(Y);
mm   = size(T,2);
a    = zeros(mm,1);
dF   = 1;
QQ   = R*Q*transpose(R);
t    = 0;
oldK = Inf;
lik  = zeros(smpl,1);
LIK  = Inf;
notsteady   = 1;

while rank(Pinf,kalman_tol) && (t<smpl)
    t  = t+1;
    d_index = data_index{t};
    if isempty(d_index)
        a = T*a;
        Pstar = T*Pstar*transpose(T)+QQ;
        Pinf  = T*Pinf*transpose(T);     
    else
        ZZ = Z(d_index,:);
        v  = Y(d_index,t)-ZZ*a;
        Finf  = ZZ*Pinf*ZZ';
        if rcond(Finf) < kalman_tol 
            if ~all(abs(Finf(:)) < kalman_tol)
                % The univariate diffuse kalman filter shoudl be used.
                return
            else
                if ~isscalar(H)                        % => Errors in the measurement equation.
                    Fstar = ZZ*Pstar*ZZ' + H(d_index,d_index); 
                else% => 
                    % case 1. No errors in the measurement (H=0) and more than one variable is observed in this state space model. 
                    % case 2. Errors in the measurement equation, but only one variable is observed in this state-space model.
                    Fstar = ZZ*Pstar*ZZ' + H;
                end
                if rcond(Fstar) < kalman_tol
                    if ~all(abs(Fstar(:))<kalman_tol)
                        % The univariate diffuse kalman filter should be used.
                        return
                    else
                        a = T*a;
                        Pstar = T*Pstar*transpose(T)+QQ;
                        Pinf  = T*Pinf*transpose(T);
                    end
                else
                    iFstar = inv(Fstar);
                    dFstar = det(Fstar);
                    Kstar  = Pstar*ZZ'*iFstar;
                    lik(t) = log(dFstar) + v'*iFstar*v + length(d_index)*log(2*pi);
                    Pinf   = T*Pinf*transpose(T);
                    Pstar  = T*(Pstar-Pstar*ZZ'*Kstar')*T'+QQ;
                    a        = T*(a+Kstar*v);
                end
            end
        else
            lik(t) = log(det(Finf));
            iFinf  = inv(Finf);
            Kinf   = Pinf*ZZ'*iFinf;
            Fstar  = ZZ*Pstar*ZZ' + H;
            Kstar  = (Pstar*ZZ'-Kinf*Fstar)*iFinf;
            Pstar  = T*(Pstar-Pstar*ZZ'*Kinf'-Pinf*ZZ'*Kstar')*T'+QQ;
            Pinf   = T*(Pinf-Pinf*ZZ'*Kinf')*T';
            a        = T*(a+Kinf*v);
        end
    end
end

if t == smpl
    warning(['There isn''t enough information to estimate the initial conditions of the nonstationary variables']);                   
    LIK = NaN;
    return
end

F_singular = 1;
while notsteady && (t<smpl)
    t = t+1;
    d_index = data_index{t};
    if isempty(d_index)
        a = T*a;
        Pstar = T*Pstar*transpose(T)+QQ;      
    else
        ZZ = Z(d_index,:);
        v = Y(d_index,t)-ZZ*a;
        if ~isscalar(H)                        % => Errors in the measurement equation.
            F = ZZ*Pstar*ZZ' + H(d_index,d_index); 
        else% => 
            % case 1. No errors in the measurement (H=0) and more than one variable is observed in this state space model. 
            % case 2. Errors in the measurement equation, but only one variable is observed in this state-space model.
            F = ZZ*Pstar*ZZ' + H;
        end
        dF = det(F);
        if rcond(F) < kalman_tol
            if ~all(abs(F(:))<kalman_tol)
                return
            else
                a     = T*a;
                Pstar = T*Pstar*T'+QQ;
            end
        else
            F_singular = 0;
            iF     = inv(F);
            lik(t) = log(dF) + v'*iF*v  + length(d_index)*log(2*pi);
            K      = Pstar*ZZ'*iF;
            a      = T*(a+K*v); 
            Pstar  = T*(Pstar-K*ZZ*Pstar)*T'+QQ;
            if t>no_more_missing_observations
                notsteady = max(abs(K(:)-oldK))>riccati_tol;
                oldK = K(:);
            end        
        end
    end
end

if F_singular == 1
    warning(['The variance of the forecast error remains singular until the end of the sample'])
    LIK = NaN;
    return
end

if t < smpl
    t0 = t+1;
    while t<smpl
        t = t+1;
        v = Y(:,t)-Z*a;
        a = T*(a+K*v);
        lik(t) = v'*iF*v;
    end
    lik(t0:smpl) = lik(t0:smpl) + log(dF) + pp*log(2*pi);
end

lik = lik/2;

LIK    = sum(lik(start:end));% Minus the log-likelihood.