function dr=set_state_space(dr,M_)
% function dr = set_state_space(dr,M_)
% finds the state vector for structural state space representation
% sets many fields of dr 
%
% INPUTS
%   dr: structure of decision rules for stochastic simulations
%  
% OUTPUTS
%   dr: structure of decision rules for stochastic simulations
%  
% ALGORITHM
%   ...
% SPECIAL REQUIREMENTS
%   none
%  

% Copyright (C) 1996-2010 Dynare Team
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

max_lead = M_.maximum_endo_lead;
max_lag = M_.maximum_endo_lag;
endo_nbr = M_.endo_nbr;
lead_lag_incidence = M_.lead_lag_incidence;
klen = max_lag + max_lead + 1;

fwrd_var = find(any(lead_lag_incidence(max_lag+2:end,:),1))';
if max_lag > 0
    pred_var = find(any(lead_lag_incidence(1,:),1))';
    both_var = intersect(pred_var,fwrd_var);
    pred_var = setdiff(pred_var,both_var);
    fwrd_var = setdiff(fwrd_var,both_var);
    stat_var = setdiff([1:endo_nbr]',union(union(pred_var,both_var),fwrd_var));  % static variables
else
    pred_var = [];
    both_var = [];
    stat_var = setdiff([1:endo_nbr]',fwrd_var);
end
nboth = length(both_var);
npred = length(pred_var);
nfwrd = length(fwrd_var);
nstatic = length(stat_var);
order_var = [ stat_var(:); pred_var(:); both_var(:); fwrd_var(:)];
inv_order_var(order_var) = (1:endo_nbr);

% building kmask for z state vector in t+1
if max_lag > 0
    kmask = [];
    if max_lead > 0 
        kmask = [cumsum(flipud(lead_lag_incidence(max_lag+2:end,order_var)),1)] ;
    end
    kmask = [kmask; flipud(cumsum(lead_lag_incidence(1,order_var),1))] ;
else
    kmask = cumsum(flipud(lead_lag_incidence(max_lag+2:klen,order_var)),1) ;
end

kmask = kmask';
kmask = kmask(:);
i_kmask = find(kmask);
nd = nnz(kmask);           % size of the state vector
kmask(i_kmask) = (1:nd);
% auxiliary equations

% elements that are both in z(t+1) and z(t)
k1 = find([kmask(1:end-M_.endo_nbr) & kmask(M_.endo_nbr+1:end)] );
kad = [];
kae = [];
if ~isempty(k1)
    kad = kmask(k1+M_.endo_nbr);
    kae = kmask(k1);
end

% composition of state vector
% col 1: variable;           col 2: lead/lag in z(t+1); 
% col 3: A cols for t+1 (D); col 4: A cols for t (E)
kstate = [ repmat([1:endo_nbr]',klen-1,1) kron([klen:-1:2]',ones(endo_nbr,1)) ...
           zeros((klen-1)*endo_nbr,2)];
kiy = flipud(lead_lag_incidence(:,order_var))';
kiy = kiy(:);
if max_lead > 0
    kstate(1:endo_nbr,3) = kiy(1:endo_nbr)-nnz(lead_lag_incidence(max_lag+1,:));  
    kstate(kstate(:,3) < 0,3) = 0;
    kstate(endo_nbr+1:end,4) = kiy(2*endo_nbr+1:end);  
else
    kstate(:,4) = kiy(endo_nbr+1:end);  
end
kstate = kstate(i_kmask,:);
   
dr.order_var = order_var;
dr.inv_order_var = inv_order_var';
dr.nstatic = nstatic;
dr.npred = npred+nboth;
dr.kstate = kstate;
dr.kad = kad;
dr.kae = kae;
dr.nboth = nboth;
dr.nfwrd = nfwrd;
% number of forward variables in the state vector
dr.nsfwrd = nfwrd+nboth;
% number of predetermined variables in the state vector
dr.nspred = npred+nboth;

dr.transition_auxiliary_variables = [];
