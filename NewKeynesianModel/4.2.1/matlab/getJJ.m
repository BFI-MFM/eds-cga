function [JJ, H, gam, gp] = getJJ(A, B, M_,oo_,options_,kronflag,indx,indexo,mf,nlags,useautocorr)

% Copyright (C) 2010-2011 Dynare Team
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

if nargin<7 | isempty(indx), indx = [1:M_.param_nbr];, end,
if nargin<8 | isempty(indexo), indexo = [];, end,
if nargin<10 | isempty(nlags), nlags=3; end,
if nargin<11 | isempty(useautocorr), useautocorr=0; end,

%   if useautocorr,
warning('off','MATLAB:divideByZero')
%   end
if kronflag == -1,
    fun = 'thet2tau';
    params0 = M_.params;
    JJ = fjaco(fun,[sqrt(diag(M_.Sigma_e(indexo,indexo))); M_.params(indx)],indx,indexo,1,mf,nlags,useautocorr);
    M_.params = params0;
    params0 = M_.params;
    H = fjaco(fun,[sqrt(diag(M_.Sigma_e(indexo,indexo))); M_.params(indx)],indx,indexo,0,mf,nlags,useautocorr);
    M_.params = params0;
    assignin('base','M_', M_);
    assignin('base','oo_', oo_);
else
    [H, dA, dOm, dYss, gp] = getH(A, B, M_,oo_,kronflag,indx,indexo);
    gp = reshape(gp,size(gp,1)*size(gp,2),size(gp,3));
    gp = [dYss; gp];
    %   if isempty(H),
    %     JJ = [];
    %     GAM = [];
    %     return
    %   end
    m = length(A);
    GAM =  lyapunov_symm(A,B*M_.Sigma_e*B',options_.qz_criterium,options_.lyapunov_complex_threshold,1);
    k = find(abs(GAM) < 1e-12);
    GAM(k) = 0;
    %   if useautocorr,
    sdy = sqrt(diag(GAM));
    sy = sdy*sdy';
    %   end
    
    %   BB = dOm*0;
    %   for j=1:length(indx),
    %     BB(:,:,j)= dA(:,:,j)*GAM*A'+A*GAM*dA(:,:,j)'+dOm(:,:,j);
    %   end
    %   XX =  lyapunov_symm_mr(A,BB,options_.qz_criterium,options_.lyapunov_complex_threshold,0);
    for j=1:length(indexo),
        dum =  lyapunov_symm(A,dOm(:,:,j),options_.qz_criterium,options_.lyapunov_complex_threshold,2);
        %     dum =  XX(:,:,j);
        k = find(abs(dum) < 1e-12);
        dum(k) = 0;
        if useautocorr
            dsy = 1/2./sdy.*diag(dum);
            dsy = dsy*sdy'+sdy*dsy';
            dum1=dum;
            dum1 = (dum1.*sy-dsy.*GAM)./(sy.*sy);
            dum1 = dum1-diag(diag(dum1))+diag(diag(dum));
            dumm = dyn_vech(dum1(mf,mf));
        else
            dumm = dyn_vech(dum(mf,mf));
        end
        for i=1:nlags,
            dum1 = A^i*dum;
            if useautocorr
                dum1 = (dum1.*sy-dsy.*(A^i*GAM))./(sy.*sy);
            end
            dumm = [dumm; vec(dum1(mf,mf))];
        end
        JJ(:,j) = dumm;
    end
    nexo = length(indexo);
    for j=1:length(indx),
        dum =  lyapunov_symm(A,dA(:,:,j+nexo)*GAM*A'+A*GAM*dA(:,:,j+nexo)'+dOm(:,:,j+nexo),options_.qz_criterium,options_.lyapunov_complex_threshold,2);
        %     dum =  XX(:,:,j);
        k = find(abs(dum) < 1e-12);
        dum(k) = 0;
        if useautocorr
            dsy = 1/2./sdy.*diag(dum);
            dsy = dsy*sdy'+sdy*dsy';
            dum1=dum;
            dum1 = (dum1.*sy-dsy.*GAM)./(sy.*sy);
            dum1 = dum1-diag(diag(dum1))+diag(diag(dum));
            dumm = dyn_vech(dum1(mf,mf));
        else
            dumm = dyn_vech(dum(mf,mf));
        end
        for i=1:nlags,
            dum1 = A^i*dum;
            for ii=1:i,
                dum1 = dum1 + A^(ii-1)*dA(:,:,j+nexo)*A^(i-ii)*GAM;
            end
            if useautocorr
                dum1 = (dum1.*sy-dsy.*(A^i*GAM))./(sy.*sy);
            end
            dumm = [dumm; vec(dum1(mf,mf))];
        end
        JJ(:,j+nexo) = dumm;
    end
    
    JJ = [ [zeros(length(mf),nexo) dYss(mf,:)]; JJ];
    
end
if nargout >2,
    %     sy=sy(mf,mf);
    options_.ar=nlags;
    [GAM,stationary_vars] = th_autocovariances(oo_.dr,oo_.dr.order_var(mf),M_,options_);
    sy=sqrt(diag(GAM{1}));
    sy=sy*sy';
    if useautocorr,
        sy=sy-diag(diag(sy))+eye(length(mf));
        GAM{1}=GAM{1}./sy;
    else
        for j=1:nlags,
            GAM{j+1}=GAM{j+1}.*sy;
        end
    end
    gam = dyn_vech(GAM{1});
    for j=1:nlags,
        gam = [gam; vec(GAM{j+1})];
    end
end
gam = [oo_.dr.ys(oo_.dr.order_var(mf)); gam];

%   if useautocorr,
warning('on','MATLAB:divideByZero')
%   end