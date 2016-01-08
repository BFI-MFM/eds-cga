function [alphahat,epsilonhat,etahat,atilde,P,aK,PK,decomp] = missing_DiffuseKalmanSmootherH1_Z(T,Z,R,Q,H,Pinf1,Pstar1,Y,pp,mm,smpl,data_index,nk,kalman_tol,decomp_flag)

% function [alphahat,epsilonhat,etahat,a,aK,PK,decomp] = DiffuseKalmanSmoother1(T,Z,R,Q,H,Pinf1,Pstar1,Y,pp,mm,smpl,data_index,nk,kalman_tol,decomp_flag)
% Computes the diffuse kalman smoother without measurement error, in the case of a non-singular var-cov matrix 
%
% INPUTS
%    T:        mm*mm matrix
%    Z:        pp*mm matrix  
%    R:        mm*rr matrix
%    Q:        rr*rr matrix
%    H:        pp*pp matrix variance of measurement errors    
%    Pinf1:    mm*mm diagonal matrix with with q ones and m-q zeros
%    Pstar1:   mm*mm variance-covariance matrix with stationary variables
%    Y:        pp*1 vector
%    pp:       number of observed variables
%    mm:       number of state variables
%    smpl:     sample size
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    nk        number of forecasting periods
%    kalman_tol   tolerance for reciprocal condition number
%    decomp_flag  if true, compute filter decomposition
%             
% OUTPUTS
%    alphahat: smoothed variables (a_{t|T})
%    epsilonhat:smoothed measurement errors
%    etahat:   smoothed shocks
%    atilde:   matrix of updated variables (a_{t|t})
%    aK:       3D array of k step ahead filtered state variables (a_{t+k|t)
%              (meaningless for periods 1:d)
%    P:        3D array of one-step ahead forecast error variance
%              matrices
%    PK:       4D array of k-step ahead forecast error variance
%              matrices (meaningless for periods 1:d)
%    decomp:   decomposition of the effect of shocks on filtered values
%  
% SPECIAL REQUIREMENTS
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series 
%   Analysis, vol. 24(1), pp. 85-98). 

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

% modified by M. Ratto:
% new output argument aK (1-step to k-step predictions)
% new options_.nk: the max step ahed prediction in aK (default is 4)
% new crit1 value for rank of Pinf
% it is assured that P is symmetric 

d = 0;
decomp = [];
spinf           = size(Pinf1);
spstar          = size(Pstar1);
v               = zeros(pp,smpl);
a               = zeros(mm,smpl+1);
atilde          = zeros(mm,smpl);
aK              = zeros(nk,mm,smpl+nk);
PK              = zeros(nk,mm,mm,smpl+nk);
iF              = zeros(pp,pp,smpl);
Fstar           = zeros(pp,pp,smpl);
iFinf           = zeros(pp,pp,smpl);
K               = zeros(mm,pp,smpl);
L               = zeros(mm,mm,smpl);
Linf            = zeros(mm,mm,smpl);
Kstar           = zeros(mm,pp,smpl);
P               = zeros(mm,mm,smpl+1);
Pstar           = zeros(spstar(1),spstar(2),smpl+1); Pstar(:,:,1) = Pstar1;
Pinf            = zeros(spinf(1),spinf(2),smpl+1); Pinf(:,:,1) = Pinf1;
crit1       = 1.e-8;
steady          = smpl;
rr              = size(Q,1);
QQ              = R*Q*transpose(R);
QRt             = Q*transpose(R);
alphahat        = zeros(mm,smpl);
etahat          = zeros(rr,smpl);
epsilonhat          = zeros(rr,smpl);
r               = zeros(mm,smpl+1);

t = 0;
while rank(Pinf(:,:,t+1),crit1) && t<smpl
    t = t+1;
    di = data_index{t};
    if isempty(di)
        atilde(:,t) = a(:,t);
        Linf(:,:,t)     = T;
        Pstar(:,:,t+1)  = T*Pstar(:,:,t)*T' + QQ;
        Pinf(:,:,t+1)   = T*Pinf(:,:,t)*T';
    else
        ZZ = Z(di,:);
        v(di,t)= Y(di,t) - ZZ*a(:,t);
        Finf = ZZ*Pinf(:,:,t)*ZZ';
        if rcond(Finf) < kalman_tol
            if ~all(abs(Finf(:)) < kalman_tol)
                % The univariate diffuse kalman filter should be used.
                return
            else
                Fstar(:,:,t)  = ZZ*Pstar(:,:,t)*ZZ' + H(di,di);
                if rcond(Fstar(:,:,t)) < kalman_tol
                    if ~all(abs(Fstar(:,:,t))<kalman_tol)
                        % The univariate diffuse kalman filter should be used.
                        return
                    else
                        a(:,:,t+1) = T*a(:,:,t);
                        Pstar(:,:,t+1) = T*Pstar(:,:,t)*transpose(T)+QQ;
                        Pinf(:,:,t+1)  = T*Pinf(:,:,t)*transpose(T);
                    end
                else
                    iFstar = inv(Fstar(:,:,t));
                    Kstar(:,:,t)  = Pstar(:,:,t)*ZZ'*iFstar(:,:,t);
                    Pinf(:,:,t+1)   = T*Pinf(:,:,t)*transpose(T);
                    Pstar(:,:,t+1)  = T*(Pstar(:,:,t)-Pstar(:,:,t)*ZZ'*Kstar(:,:,t)')*T'+QQ;
                    a(:,:,t+1)        = T*(a(:,:,t)+Kstar(:,:,t)*v(:,t));
                end
            end
        else
            iFinf(di,di,t)  = inv(Finf);
            PZI             = Pinf(:,:,t)*ZZ'*iFinf(di,di,t);
            atilde(:,t)     = a(:,t) + PZI*v(di,t);
            Kinf(:,di,t)    = T*PZI;
            Linf(:,:,t)     = T - Kinf(:,di,t)*ZZ;
            Fstar(di,di,t)  = ZZ*Pstar(:,:,t)*ZZ' + H(di,di);
            Kstar(:,di,t)   = (T*Pstar(:,:,t)*ZZ'-Kinf(:,di,t)*Fstar(di,di,t))*iFinf(di,di,t);
            Pstar(:,:,t+1)  = T*Pstar(:,:,t)*T'-T*Pstar(:,:,t)*ZZ'*Kinf(:,di,t)'-T*Pinf(:,:,t)*ZZ'*Kstar(:,di,t)' + QQ;
            Pinf(:,:,t+1)   = T*Pinf(:,:,t)*T'-T*Pinf(:,:,t)*ZZ'*Kinf(:,di,t)';
        end
        a(:,t+1)    = T*atilde(:,t);
        aK(1,:,t+1)         = a(:,t+1);
        % isn't a meaningless as long as we are in the diffuse part? MJ
        for jnk=2:nk,
            aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
        end
    end
end
d = t;
P(:,:,d+1) = Pstar(:,:,d+1);
iFinf = iFinf(:,:,1:d);
Linf  = Linf(:,:,1:d);
Fstar = Fstar(:,:,1:d);
Kstar = Kstar(:,:,1:d);
Pstar = Pstar(:,:,1:d);
Pinf  = Pinf(:,:,1:d);
notsteady = 1;
while notsteady && t<smpl
    t = t+1;
    P(:,:,t)=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
    di = data_index{t};
    if isempty(di)
        atilde(:,t) = a(:,t);
        L(:,:,t)        = T;
        P(:,:,t+1)      = T*P(:,:,t)*T' + QQ;
    else
        ZZ = Z(di,:);
        v(di,t)      = Y(di,t) - ZZ*a(:,t);
        F = ZZ*P(:,:,t)*ZZ' + H(di,di);
        if rcond(F) < kalman_tol
            return              
        end    
        iF(di,di,t)   = inv(F);
        PZI         = P(:,:,t)*ZZ'*iF(di,di,t);
        atilde(:,t) = a(:,t) + PZI*v(di,t);
        K(:,di,t)    = T*PZI;
        L(:,:,t)    = T-K(:,di,t)*ZZ;
        P(:,:,t+1)  = T*P(:,:,t)*T'-T*P(:,:,t)*ZZ'*K(:,di,t)' + QQ;
    end
    a(:,t+1)    = T*atilde(:,t);
    Pf          = P(:,:,t);
    aK(1,:,t+1) = a(:,t+1);
    for jnk=1:nk
        Pf = T*Pf*T' + QQ;
        PK(jnk,:,:,t+jnk) = Pf;
        if jnk>1
            aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
        end
    end
    %    notsteady   = ~(max(max(abs(P(:,:,t+1)-P(:,:,t))))<kalman_tol);
end
% $$$ if t<smpl
% $$$     PZI_s = PZI;
% $$$     K_s = K(:,:,t);
% $$$     iF_s = iF(:,:,t);
% $$$     P_s = P(:,:,t+1);
% $$$     P  = cat(3,P(:,:,1:t),repmat(P_s,[1 1 smpl-t]));
% $$$     iF = cat(3,iF(:,:,1:t),repmat(iF_s,[1 1 smpl-t]));
% $$$     L  = cat(3,L(:,:,1:t),repmat(T-K_s*Z,[1 1 smpl-t]));
% $$$     K  = cat(3,K(:,:,1:t),repmat(T*P_s*Z'*iF_s,[1 1 smpl-t]));
% $$$ end
% $$$ while t<smpl
% $$$     t=t+1;
% $$$     v(:,t) = Y(:,t) - Z*a(:,t);
% $$$     atilde(:,t) = a(:,t) + PZI*v(:,t);
% $$$     a(:,t+1) = T*atilde(:,t);
% $$$     Pf          = P(:,:,t);
% $$$     for jnk=1:nk,
% $$$   Pf = T*Pf*T' + QQ;
% $$$         aK(jnk,:,t+jnk) = T^jnk*atilde(:,t);
% $$$   PK(jnk,:,:,t+jnk) = Pf;
% $$$     end
% $$$ end
t = smpl+1;
while t>d+1
    t = t-1;
    di = data_index{t};
    if isempty(di)
        r(:,t) = L(:,:,t)'*r(:,t+1);
    else
        ZZ = Z(di,:);
        r(:,t) = ZZ'*iF(di,di,t)*v(di,t) + L(:,:,t)'*r(:,t+1);
    end
    alphahat(:,t)       = a(:,t) + P(:,:,t)*r(:,t);
    etahat(:,t) = QRt*r(:,t);
end
if d
    r0 = zeros(mm,d+1); 
    r0(:,d+1) = r(:,d+1);
    r1 = zeros(mm,d+1);
    for t = d:-1:1
        r0(:,t) = Linf(:,:,t)'*r0(:,t+1);
        di = data_index{t};
        if isempty(di)
            r1(:,t) = Linf(:,:,t)'*r1(:,t+1);
        else
            r1(:,t) = Z(di,:)'*(iFinf(di,di,t)*v(di,t)-Kstar(:,di,t)'*r0(:,t+1)) ...
                      + Linf(:,:,t)'*r1(:,t+1);
        end
        alphahat(:,t)   = a(:,t) + Pstar(:,:,t)*r0(:,t) + Pinf(:,:,t)*r1(:,t);
        etahat(:,t)             = QRt*r0(:,t);
    end
end

if decomp_flag
    decomp = zeros(nk,mm,rr,smpl+nk);
    ZRQinv = inv(Z*QQ*Z');
    for t = max(d,1):smpl
        di = data_index{t};
        % calculate eta_tm1t
        eta_tm1t = QRt*Z(di,:)'*iF(di,di,t)*v(di,t);
        AAA = P(:,:,t)*Z(di,:)'*ZRQinv(di,di)*bsxfun(@times,Z(di,:)*R,eta_tm1t');
        % calculate decomposition
        Ttok = eye(mm,mm); 
        decomp(1,:,:,t+1) = AAA;
        for h = 2:nk
            AAA = T*AAA;
            decomp(h,:,:,t+h) = AAA;
        end
    end
end

epsilonhat = Y-Z*alphahat;
