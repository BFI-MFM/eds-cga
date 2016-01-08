function [J,M_,dr] = dyn_ramsey_dynamic_(ys,lbar,M_,options_,dr,it_)
% function J = dyn_ramsey_dynamic_(ys,lbar)
% dyn_ramsey_dynamic_ sets up the Jacobian of the expanded model for optimal
% policies. It modifies several fields of M_
%
% INPUTS:
%     ys:         steady state of original endogenous variables
%     lbar:       steady state of Lagrange multipliers
%
% OUPTUTS: 
%     J:          jaocobian of expanded model
%  
% SPECIAL REQUIREMENTS
%     none

% Copyright (C) 2003-2011 Dynare Team
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

persistent old_lead_lag

% retrieving model parameters
endo_nbr = M_.endo_nbr;
i_endo_nbr = 1:endo_nbr;
endo_names = M_.endo_names;
%  exo_nbr = M_.exo_nbr+M_.exo_det_nbr;
%  exo_names = vertcat(M_.exo_names,M_.exo_det_names);
exo_nbr = M_.exo_nbr;
exo_names = M_.exo_names;
i_leadlag = M_.lead_lag_incidence;
max_lead = M_.maximum_lead;
max_endo_lead = M_.maximum_endo_lead;
max_lag = M_.maximum_lag;
max_endo_lag = M_.maximum_endo_lag;
leadlag_nbr = max_lead+max_lag+1;
fname = M_.fname;
% instr_names = options_.olr_inst;
% instr_nbr =  size(options_.olr_inst,1);

% discount factor
beta = options_.planner_discount;

% storing original values
orig_model.endo_nbr = endo_nbr;
orig_model.orig_endo_nbr = M_.orig_endo_nbr;
orig_model.aux_vars = M_.aux_vars;
orig_model.endo_names = endo_names;
orig_model.lead_lag_incidence = i_leadlag;
orig_model.maximum_lead = max_lead;
orig_model.maximum_endo_lead = max_endo_lead;
orig_model.maximum_lag = max_lag;
orig_model.maximum_endo_lag = max_endo_lag;

y = repmat(ys,1,max_lag+max_lead+1);
k = find(i_leadlag');

% retrieving derivatives of the objective function
[U,Uy,Uyy] = feval([fname '_objective_static'],ys,zeros(1,exo_nbr), M_.params);
Uy = Uy';
Uyy = reshape(Uyy,endo_nbr,endo_nbr);

% retrieving derivatives of original model
[f,fJ,fH] = feval([fname '_dynamic'],y(k),zeros(2,exo_nbr), M_.params, ys, ...
                  it_);
instr_nbr = endo_nbr - size(f,1);
mult_nbr = endo_nbr-instr_nbr;

% parameters for expanded model
endo_nbr1 = 2*endo_nbr-instr_nbr+exo_nbr;
max_lead1 = max_lead + max_lag;
max_lag1 = max_lead1;
max_leadlag1 = max_lead1;

% adding new variables names
endo_names1 = endo_names;
% adding shocks to endogenous variables
endo_names1 = char(endo_names1, exo_names);
% adding multipliers names
for i=1:mult_nbr;
    endo_names1 = char(endo_names1,['mult_' int2str(i)]);
end

% expanding matrix of lead/lag incidence
%
% multipliers lead/lag incidence
i_mult = [];
for i=1:leadlag_nbr
    i_mult = [any(fJ(:,nonzeros(i_leadlag(i,:))) ~= 0,2)' ; i_mult];
end
% putting it all together:
% original variables, exogenous variables made endogenous, multipliers
%
% number of original dynamic variables 
n_dyn = nnz(i_leadlag);
% numbering columns of dynamic multipliers to be put in the last columns
% of the new Jacobian
i_leadlag1 = [cumsum(i_leadlag(1:max_lag,:),1); ...
              repmat(i_leadlag(max_lag+1,:),leadlag_nbr,1); ...
              flipud(cumsum(flipud(i_leadlag(max_lag+2:end,:)),1))];
i_leadlag1 = i_leadlag1';
k = find(i_leadlag1 > 0);
n = length(k);
i_leadlag1(k) = 1:n;
i_leadlag1 = i_leadlag1';
i_mult = i_mult';
k = find(i_mult > 0);
i_mult(k) = n+leadlag_nbr*exo_nbr+(1:length(k));
i_mult = i_mult';
i_leadlag1 = [  i_leadlag1 ...
                [zeros(max_lag,exo_nbr);...
                 reshape(n+(1:leadlag_nbr*exo_nbr),exo_nbr,leadlag_nbr)'; ...
                 zeros(max_lead,exo_nbr)] ...
                [zeros(max_lag,mult_nbr);...
                 i_mult;...
                 zeros(max_lead,mult_nbr)]];
i_leadlag1 = i_leadlag1';
k = find(i_leadlag1 > 0);
n = length(k);
i_leadlag1(k) = 1:n;
i_leadlag1 = i_leadlag1';

% building Jacobian of expanded model
jacobian = zeros(endo_nbr+mult_nbr,nnz(i_leadlag1)+exo_nbr);
% derivatives of f.o.c. w.r. to endogenous variables
% to be rearranged further down
lbarfH = lbar'*fH; 
% indices of Hessian columns
n1 = nnz(i_leadlag)+exo_nbr;
iH = reshape(1:n1^2,n1,n1);
J = zeros(endo_nbr1,nnz(i_leadlag1)+exo_nbr);
% second order derivatives of objective function
J(1:endo_nbr,i_leadlag1(max_leadlag1+1,1:endo_nbr)) = Uyy;
% loop on lead/lags in expanded model 
for i=1:2*max_leadlag1 + 1
    % index of variables at the current lag in expanded model
    kc = find(i_leadlag1(i,i_endo_nbr) > 0);
    t1 = max(1,i-max_leadlag1);
    t2 = min(i,max_leadlag1+1);
    % loop on lead/lag blocks of relevant 1st order derivatives
    for j = t1:t2
        % derivatives w.r. endogenous variables
        ic =  find(i_leadlag(i-j+1,:) > 0 );
        kc1 =  i_leadlag(i-j+1,ic);
        [junk,ic1,ic2] = intersect(ic,kc);
        kc2 = i_leadlag1(i,kc(ic2));
        ir = find(i_leadlag(max_leadlag1+2-j,:) > 0 );
        kr1 = i_leadlag(max_leadlag1+2-j,ir);
        J(ir,kc2) = J(ir,kc2) + beta^(j-max_lead-1)...
            *reshape(lbarfH(iH(kr1,kc1)),length(kr1),length(kc1));
    end
end
% derivatives w.r. aux. variables for lead/lag exogenous shocks
for i=1:leadlag_nbr
    kc = i_leadlag1(max_lag+i,endo_nbr+(1:exo_nbr));
    ir = find(i_leadlag(leadlag_nbr+1-i,:) > 0);
    kr1 = i_leadlag(leadlag_nbr+1-i,ir);
    J(ir,kc) = beta^(i-max_lead-1)...
        *reshape(lbarfH(iH(kr1,n_dyn+(1:exo_nbr))),length(kr1), ...
                 exo_nbr);
end
% derivatives w.r. Lagrange multipliers
for i=1:leadlag_nbr
    ic1 = find(i_leadlag(leadlag_nbr+1-i,:) > 0);
    kc1 = i_leadlag(leadlag_nbr+1-i,ic1);
    ic2 = find(i_leadlag1(max_lag+i,endo_nbr+exo_nbr+(1:mult_nbr)) > 0);
    kc2 = i_leadlag1(max_lag+i,endo_nbr+exo_nbr+ic2);
    J(ic1,kc2) = beta^(i-max_lead-1)*fJ(ic2,kc1)';
end

% Jacobian of original equations
%
% w.r. endogenous variables
ir = endo_nbr+(1:endo_nbr-instr_nbr);
for i=1:leadlag_nbr
    ic1 = find(i_leadlag(i,:) > 0);
    kc1 = i_leadlag(i,ic1);
    ic2 = find(i_leadlag1(max_lead+i,:) > 0);
    kc2 = i_leadlag1(max_lead+i,ic2);
    [junk,junk,ic3] = intersect(ic1,ic2);
    J(ir,kc2(ic3)) = fJ(:,kc1);
end
% w.r. exogenous variables
J(ir,nnz(i_leadlag1)+(1:exo_nbr)) = fJ(:,nnz(i_leadlag)+(1:exo_nbr));

% auxiliary variable for exogenous shocks
ir = 2*endo_nbr-instr_nbr+(1:exo_nbr);
kc = i_leadlag1(leadlag_nbr,endo_nbr+(1:exo_nbr));
J(ir,kc) = eye(exo_nbr);
J(ir,nnz(i_leadlag1)+(1:exo_nbr)) = -eye(exo_nbr);

% eliminating empty columns

% getting indices of nonzero entries
m = find(i_leadlag1');
n1 = max_lag1*endo_nbr1+1;
n2 = n1+endo_nbr-1;


n = length(m);
k = 1:size(J,2);

for i=1:n
    if sum(abs(J(:,i))) < 1e-8
        if m(i) < n1 || m(i) > n2
            k(i) = 0;
            m(i) = 0;
        end
    end
end

J = J(:,nonzeros(k));
i_leadlag1 = zeros(size(i_leadlag1))';
i_leadlag1(nonzeros(m)) = 1:nnz(m);
i_leadlag1 = i_leadlag1';

%eliminating lags in t-2 and leads in t+2, if possible
if all(i_leadlag1(5,:)==0)
    i_leadlag1 = i_leadlag1(1:4,:);
    max_lead1 = 1;
end

if all(i_leadlag1(1,:)==0)
    i_leadlag1 = i_leadlag1(2:4,:);
    max_lag1 = 1;
end

% setting expanded model parameters
% storing original values
M_.endo_nbr = endo_nbr1;
% Consider that there is no auxiliary variable, because otherwise it
% interacts badly with the auxiliary variables from the preprocessor.
M_.orig_endo_nbr = endo_nbr1;
M_.aux_vars = [];
M_.endo_names = endo_names1;
M_.lead_lag_incidence = i_leadlag1;
M_.maximum_lead = max_lead1;
M_.maximum_endo_lead = max_lead1;
M_.maximum_lag = max_lag1;
M_.maximum_endo_lag = max_lag1;
M_.orig_model = orig_model;

if isfield(options_,'varobs') && (any(size(i_leadlag1,2) ~= size(old_lead_lag,2)) || any(any(i_leadlag1 ~= old_lead_lag)))
    global bayestopt_
    dr = set_state_space(dr,M_);
    nstatic = dr.nstatic;          % Number of static variables. 
    npred = dr.npred;              % Number of predetermined variables.

    var_obs_index = [];
    k1 = [];
    for i=1:size(options_.varobs,1);
        var_obs_index = [var_obs_index strmatch(deblank(options_.varobs(i,:)),M_.endo_names(dr.order_var,:),'exact')];
        k1 = [k1 strmatch(deblank(options_.varobs(i,:)),M_.endo_names, 'exact')];
    end
    % Define union of observed and state variables
    k2 = union(var_obs_index',[nstatic+1:nstatic+npred]');
    % Set restrict_state to postion of observed + state variables in expanded state vector.
    dr.restrict_var_list = k2;
    [junk,ic] = intersect(k2,nstatic+(1:npred)');
    dr.restrict_columns = ic;
    % set mf0 to positions of state variables in restricted state vector for likelihood computation.
    [junk,bayestopt_.mf0] = ismember([dr.nstatic+1:dr.nstatic+dr.npred]',k2);
    % Set mf1 to positions of observed variables in restricted state vector for likelihood computation.
    [junk,bayestopt_.mf1] = ismember(var_obs_index,k2); 
    % Set mf2 to positions of observed variables in expanded state vector for filtering and smoothing.
    bayestopt_.mf2  = var_obs_index;
    bayestopt_.mfys = k1;
    old_lead_lag = i_leadlag1;
end