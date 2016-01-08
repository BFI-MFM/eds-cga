function [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK,T,R,P,PK,decomp] = DsgeSmoother(xparam1,gend,Y,data_index,missing_value)
% Estimation of the smoothed variables and innovations. 
% 
% INPUTS 
%   o xparam1       [double]   (p*1) vector of (estimated) parameters. 
%   o gend          [integer]  scalar specifying the number of observations ==> varargin{1}.
%   o data          [double]   (T*n) matrix of data.
%   o data_index    [cell]      1*smpl cell of column vectors of indices.
%   o missing_value 1 if missing values, 0 otherwise
%  
% OUTPUTS 
%   o alphahat      [double]  (m*T) matrix, smoothed endogenous variables.
%   o etahat        [double]  (r*T) matrix, smoothed structural shocks (r>n is the umber of shocks).
%   o epsilonhat    [double]  (n*T) matrix, smoothed measurement errors.
%   o ahat          [double]  (m*T) matrix, one step ahead filtered (endogenous) variables.
%   o SteadyState   [double]  (m*1) vector specifying the steady state level of each endogenous variable.
%   o trend_coeff   [double]  (n*1) vector, parameters specifying the slope of the trend associated to each observed variable.
%   o aK            [double]  (K,n,T+K) array, k (k=1,...,K) steps ahead filtered (endogenous) variables.
%   o T and R       [double]  Matrices defining the state equation (T is the (m*m) transition matrix).
%    P:             3D array of one-step ahead forecast error variance
%                   matrices
%    PK:            4D array of k-step ahead forecast error variance
%                   matrices (meaningless for periods 1:d)
%    
% ALGORITHM 
%   Diffuse Kalman filter (Durbin and Koopman)       
%
% SPECIAL REQUIREMENTS
%   None

% Copyright (C) 2006-2011 Dynare Team
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

global bayestopt_ M_ oo_ estim_params_ options_

alphahat        = [];
etahat  = [];
epsilonhat      = [];
ahat          = [];
SteadyState   = [];
trend_coeff   = [];
aK            = [];
T             = [];
R             = [];
P             = [];
PK            = [];
decomp        = [];
nobs            = size(options_.varobs,1);
smpl          = size(Y,2);

set_all_parameters(xparam1);

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------
oo_.dr.restrict_var_list = bayestopt_.smoother_var_list;
oo_.dr.restrict_columns = bayestopt_.smoother_restrict_columns;
[T,R,SteadyState] = dynare_resolve('restrict');

bayestopt_.mf = bayestopt_.smoother_mf;
if options_.noconstant
    constant = zeros(nobs,1);
else
    if options_.loglinear == 1
        constant = log(SteadyState(bayestopt_.mfys));
    else
        constant = SteadyState(bayestopt_.mfys);
    end
end
trend_coeff = zeros(nobs,1);
if bayestopt_.with_trend == 1
    trend_coeff = zeros(nobs,1);
    t = options_.trend_coeffs;
    for i=1:length(t)
        if ~isempty(t{i})
            trend_coeff(i) = evalin('base',t{i});
        end
    end
    trend = constant*ones(1,gend)+trend_coeff*(1:gend);
else
    trend = constant*ones(1,gend);
end
start = options_.presample+1;
np    = size(T,1);
mf    = bayestopt_.smoother_mf;
% ------------------------------------------------------------------------------
%  3. Initial condition of the Kalman filter
% ------------------------------------------------------------------------------
% 
%  C'est ici qu'il faut déterminer Pinf et Pstar. Si le modèle est stationnaire,
%  alors il suffit de poser Pstar comme la solution de l'éuation de Lyapounov et
%  Pinf=[].
%
Q = M_.Sigma_e;
H = M_.H;

kalman_algo = options_.kalman_algo;
if options_.lik_init == 1               % Kalman filter
    if kalman_algo ~= 2
        kalman_algo = 1;
    end
    Pstar = lyapunov_symm(T,R*Q*transpose(R),options_.qz_criterium,options_.lyapunov_complex_threshold);
    Pinf        = [];
elseif options_.lik_init == 2 % Old Diffuse Kalman filter
    if kalman_algo ~= 2
        kalman_algo = 1;
    end
    Pstar = options_.Harvey_scale_factor*eye(np);
    Pinf        = [];
elseif options_.lik_init == 3 % Diffuse Kalman filter
    if kalman_algo ~= 4
        kalman_algo = 3;
    end
    [Z,ST,R1,QT,Pstar,Pinf] = schur_statespace_transformation(mf,T,R,Q,options_.qz_criterium);
end
kalman_tol = options_.kalman_tol;
riccati_tol = options_.riccati_tol;
data1 = Y-trend;
% -----------------------------------------------------------------------------
%  4. Kalman smoother
% -----------------------------------------------------------------------------
if isequal(H,0)
    H = zeros(nobs,nobs);
end

if ~missing_value
    for i=1:smpl
        data_index{i}=(1:nobs)';
    end
end

if kalman_algo == 1 || kalman_algo == 2
    ST = T;
    R1 = R;
    Z = zeros(nobs,size(T,2));
    for i=1:nobs
        Z(i,mf(i)) = 1;
    end
end

if kalman_algo == 1 || kalman_algo == 3
    [alphahat,epsilonhat,etahat,ahat,P,aK,PK,decomp] = missing_DiffuseKalmanSmootherH1_Z(ST, ...
                                                      Z,R1,Q,H,Pinf,Pstar, ...
                                                      data1,nobs,np,smpl,data_index, ...
                                                      options_.nk,kalman_tol,options_.filter_decomposition);
    if isequal(alphahat,0)
        if kalman_algo == 1
            kalman_algo = 2;
        elseif kalman_algo == 3
            kalman_algo = 4;
        else
            error('This case shouldn''t happen')
        end
    end
end

if kalman_algo == 2 || kalman_algo == 4
    if estim_params_.ncn
        ST = [ zeros(nobs,nobs) Z; zeros(np,nobs) T];
        ns = size(Q,1);
        R1 = [ eye(nobs) zeros(nobs, ns); zeros(np,nobs) R];
        Q = [H zeros(nobs,ns); zeros(ns,nobs) Q]; 
        Z = [eye(nobs) zeros(nobs, np)];
        if kalman_algo == 4
            [Z,ST,R1,QT,Pstar,Pinf] = schur_statespace_transformation((1:nobs)',ST,R1,Q,options_.qz_criterium);
        end
        
    end
    [alphahat,epsilonhat,etahat,ahat,P,aK,PK,decomp] = missing_DiffuseKalmanSmootherH3_Z(ST, ...
                                                      Z,R1,Q,diag(H), ...
                                                      Pinf,Pstar,data1,nobs,np,smpl,data_index, ...
                                                      options_.nk,kalman_tol,...
                                                      options_.filter_decomposition);
end

if kalman_algo == 3 || kalman_algo == 4
    alphahat = QT*alphahat;
    ahat = QT*ahat;
    nk = options_.nk;
    for jnk=1:nk
        aK(jnk,:,:) = QT*dynare_squeeze(aK(jnk,:,:));
        for i=1:size(PK,4)
            PK(jnk,:,:,i) = QT*dynare_squeeze(PK(jnk,:,:,i))*QT';
        end
        if options_.filter_decomposition
            for i=1:size(decomp,4)
                decomp(jnk,:,:,i) = QT*dynare_squeeze(decomp(jnk,:,:,i));
            end
        end
    end
    for i=1:size(P,4)
        P(:,:,i) = QT*dynare_squeeze(P(:,:,i))*QT';
    end
end

if estim_params_.ncn && (kalman_algo == 2 || kalman_algo == 4)
    % extracting measurement errors
    % removing observed variables from the state vector
    k = nobs+(1:np);
    alphahat = alphahat(k,:);
    ahat = ahat(k,:);
    aK = aK(:,k,:,:);
    if ~isempty(PK)
        PK = PK(:,k,k,:);
    end
    if ~isempty(decomp)
        decomp = decomp(:,k,:,:);
    end
    if ~isempty(P)
        P = P(k,k,:);
    end
end
