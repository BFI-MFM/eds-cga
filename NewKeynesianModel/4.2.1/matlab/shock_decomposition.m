function oo_ = shock_decomposition(M_,oo_,options_,varlist)
% function z = shock_decomposition(M_,oo_,options_,varlist)
% Computes shocks contribution to a simulated trajectory
%
% INPUTS
%    M_:          [structure]  Definition of the model
%    oo_:         [structure]  Storage of results
%    options_:    [structure]  Options
%    varlist:     [char]       List of variables
%
% OUTPUTS
%    oo_:         [structure]  Storage of results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2009-2011 Dynare Team
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

% indices of endogenous variables
if size(varlist,1) == 0
    varlist = M_.endo_names(1:M_.orig_endo_nbr,:);
end

[i_var,nvar] = varlist_indices(varlist,M_.endo_names);

% number of variables
endo_nbr = M_.endo_nbr;

% number of shocks
nshocks = M_.exo_nbr;

% parameter set
parameter_set = options_.parameter_set;
if isempty(parameter_set)
    if isfield(oo_,'posterior_mean')
        parameter_set = 'posterior_mean';
    elseif isfield(oo_,'posterior_mode')
        parameter_set = 'posterior_mode';
    else
        error(['shock_decomposition: option parameter_set is not specified ' ...
               'and posterior mode is not available'])
    end
end

oo = evaluate_smoother(parameter_set);

% reduced form
dr = oo.dr;

% data reordering
order_var = dr.order_var;
inv_order_var = dr.inv_order_var;


% coefficients
A = dr.ghx;
B = dr.ghu;

% initialization
for i=1:nshocks
    epsilon(i,:) = eval(['oo.SmoothedShocks.' M_.exo_names(i,:)]);
end
gend = size(epsilon,2);

z = zeros(endo_nbr,nshocks+2,gend);
for i=1:endo_nbr
    z(i,end,:) = eval(['oo.SmoothedVariables.' M_.endo_names(i,:)]);
end

maximum_lag = M_.maximum_lag;
lead_lag_incidence = M_.lead_lag_incidence;

k2 = dr.kstate(find(dr.kstate(:,2) <= maximum_lag+1),[1 2]);
i_state = order_var(k2(:,1))+(min(i,maximum_lag)+1-k2(:,2))*M_.endo_nbr;
for i=1:gend
    if i > 1 & i <= maximum_lag+1
        lags = min(i-1,maximum_lag):-1:1;
    end
    
    if i > 1
        tempx = permute(z(:,1:nshocks,lags),[1 3 2]);
        m = min(i-1,maximum_lag);
        tempx = [reshape(tempx,endo_nbr*m,nshocks); zeros(endo_nbr*(maximum_lag-i+1),nshocks)];
        z(:,1:nshocks,i) = A(inv_order_var,:)*tempx(i_state,:);
        lags = lags+1;
    end

    z(:,1:nshocks,i) = z(:,1:nshocks,i) + B(inv_order_var,:).*repmat(epsilon(:,i)',endo_nbr,1);
    z(:,nshocks+1,i) = z(:,nshocks+2,i) - sum(z(:,1:nshocks,i),2);
end


oo_.shock_decomposition = z;

options_.initial_date.freq = 1;
options_.initial_date.period = 1;
options_.initial_date.sub_period = 0;

graph_decomp(z,M_.exo_names,M_.endo_names,i_var,options_.initial_date)
