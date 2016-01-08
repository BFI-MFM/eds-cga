function oo1 = evaluate_planner_objective(dr,M,oo,options)

%function oo1 = evaluate_planner_objective(dr,M,oo,options)
%  computes value of planner objective function     
% 
% INPUTS
%   dr:       (structure) decision rule
%   M:        (structure) model description
%   oo:       (structure) output results
%   options:  (structure) options
%    
% OUTPUTS
%   oo1:      (structure) updated output results
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2007-2011 Dynare Team
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

oo1 = oo;
endo_nbr = M.endo_nbr;
exo_nbr = M.exo_nbr;
nstatic = dr.nstatic;
npred = dr.npred;
lead_lag_incidence = M.lead_lag_incidence;
if (size(lead_lag_incidence,1) == 5) ...
        && (nnz(lead_lag_incidence(1,:)) > 0 ...
            || nnz(lead_lag_incidence(5,:)) > 0)
    error(['procedure not yet implemented for leads and lags on more ' ...
           'than one period'])
end
if options.ramsey_policy
    i_org = (1:M.orig_model.endo_nbr)';
else
    i_org = (1:M.endo_nbr)';
end
ipred = find(lead_lag_incidence(M.maximum_lag,:))';
orig_ifrwrd = find(M.orig_model.lead_lag_incidence(3,:));
orig_npred = nnz(M.orig_model.lead_lag_incidence(1,:));
order_var = dr.order_var;
ys = oo.dr.ys;
[U,Uy,Uyy] = feval([M.fname '_objective_static'],ys,zeros(1,exo_nbr), ...
                   M.params);
z = repmat(ys(1:M.orig_model.endo_nbr),1,3);
z = z(find(M.orig_model.lead_lag_incidence'));
[F,Fy,Fyy] = feval([M.fname '_dynamic'],z,zeros(3,exo_nbr), ...
                   M.params,2);
mu_ss = oo.dr.ys(M.orig_model.endo_nbr+exo_nbr+(1:size(F,1)));
%K = reshape(1:endo_nbr^2,endo_nbr,endo_nbr);
%K = K(order_var,order_var);
%Uyy = Uyy(:,K(:));
beta = options.planner_discount;
Gy = dr.ghx(nstatic+(1:npred),:);
Gu = dr.ghu(nstatic+(1:npred),:);
gy(dr.order_var,:) = dr.ghx;
gu(dr.order_var,:) = dr.ghu;
gy = gy(i_org,:);
gu = gu(i_org,:);
muFy = mu_ss'*Fy;
muFyy = mu_ss'*Fyy;
Zy = [eye(orig_npred,npred); gy; gy(orig_ifrwrd,:)*Gy; zeros(exo_nbr, ...
                                                  npred)];
Zu = [zeros(orig_npred,exo_nbr); gu; gy(orig_ifrwrd,:)*Gu; eye(exo_nbr, ...
                                                  exo_nbr)];
Zyy = kron(Zy,Zy);
Zuu = kron(Zu,Zu);
Zyu = kron(Zy,Zu);

Wbar = U/(1-beta);
Wy = (Uy*gy+muFy*Zy)/(eye(npred)-beta*Gy);
Wu = Uy*gu + muFy*Zu + beta*Wy*Gu;
%    N1 = Uyy*kron(gy,gy)
%    N2 = muFyy*Zyy
format long
Wyy = (Uyy*kron(gy,gy)+muFyy*Zyy)/(eye(npred^2)-beta*kron(Gy,Gy))
Wyu = Uyy*kron(gy,gu) + muFyy*Zyu + beta*Wyy*kron(Gy,Gu);
Wuu = Uyy*kron(gu,gu) + muFyy*Zuu + beta*Wyy*kron(Gu,Gu);
Wss = beta*Wuu*vec(M.Sigma_e)/(1-beta);

if options.ramsey_policy
    yhat = [oo.endo_simul; ...
            zeros(M.endo_nbr-size(oo.endo_simul,1),1)];
else
    yhat = oo.endo_simul;
end
yhat = yhat(dr.order_var(nstatic+(1:npred)),1)-dr.ys(dr.order_var(nstatic+(1:npred)))
u = oo.exo_simul(1,:)';

planner_objective_value = Wbar+Wy*yhat+Wu*u+Wyu*kron(yhat,u) ...
    + 0.5*(Wyy*kron(yhat,yhat) + Wuu*kron(u,u)+Wss);
disp(Wyy*kron(yhat,yhat))
oo1.planner_objective_value = planner_objective_value;
if ~options.noprint
    disp(' ')
    if options.ramsey_policy
        disp(sprintf('Value of planner objective function under Ramsey policy: %16.8f', ...
                     planner_objective_value))
    else
        disp(['Value of planner objective function: ' ...
              num2str(planner_objective_value)])
        disp(' ')
    end
end