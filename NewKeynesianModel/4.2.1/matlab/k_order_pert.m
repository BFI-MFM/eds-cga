function [dr,info] = k_order_pert(dr,M,options,oo)
% Compute decision rules using the k-order DLL from Dynare++

% Copyright (C) 2009-2010 Dynare Team
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

info = 0;

M.var_order_endo_names = M.endo_names(dr.order_var,:);

order = options.order;

switch(order)
  case 1
    [err, g_1] = k_order_perturbation(dr,M,options, ['.' ...
                        mexext]);
    mexErrCheck('k_order_perturbation', err);
    dr.g_1 = g_1;
  case 2
    [err, g_0, g_1, g_2] = k_order_perturbation(dr,M,options, ['.' ...
                        mexext]);
    mexErrCheck('k_order_perturbation', err);
    dr.g_0 = g_0;
    dr.g_1 = g_1;
    dr.g_2 = g_2;
  case 3
    [err, g_0, g_1, g_2, g_3] = k_order_perturbation(dr,M,options, ['.' ...
                        mexext]);
    mexErrCheck('k_order_perturbation', err);
    dr.g_0 = g_0;
    dr.g_1 = g_1;
    dr.g_2 = g_2;
    dr.g_3 = g_3;
  otherwise
    error('order > 3 isn''t implemented')
end

npred = dr.npred;

dr.ghx = g_1(:,1:npred);
dr.ghu = g_1(:,npred+1:end);

if options.loglinear == 1
    k = find(dr.kstate(:,2) <= M.maximum_endo_lag+1);
    klag = dr.kstate(k,[1 2]);
    k1 = dr.order_var;
    
    dr.ghx = repmat(1./dr.ys(k1),1,size(dr.ghx,2)).*dr.ghx.* ...
             repmat(dr.ys(k1(klag(:,1)))',size(dr.ghx,1),1);
    dr.ghu = repmat(1./dr.ys(k1),1,size(dr.ghu,2)).*dr.ghu;
end

if order > 1
    dr.ghs2 = 2*g_0;
    endo_nbr = M.endo_nbr;
    exo_nbr = M.exo_nbr;
    s0 = 0;
    s1 = 0;
    ghxx=zeros(endo_nbr, npred^2);
    ghxu=zeros(endo_nbr, npred*exo_nbr);
    ghuu=zeros(endo_nbr, exo_nbr^2);
    for i=1:size(g_2,2)
        if s0 < npred && s1 < npred
            ghxx(:,s0*npred+s1+1) = 2*g_2(:,i);
            if s1 > s0
                ghxx(:,s1*npred+s0+1) = 2*g_2(:,i);
            end
        elseif s0 < npred && s1 < npred+exo_nbr 
            ghxu(:,(s0*exo_nbr+s1-npred+1)) = 2*g_2(:,i);
        elseif s0 < npred+exo_nbr && s1 < npred+exo_nbr
            ghuu(:,(s0-npred)*exo_nbr+s1-npred +1) = 2*g_2(:,i);
            if s1 > s0
                ghuu(:,(s1-npred)*exo_nbr+s0-npred+1) = 2*g_2(:,i);
            end
        else
            error('dr1:k_order_perturbation:g_2','Unaccounted columns in g_2');
        end
        s1 = s1+1;
        if s1 == npred+exo_nbr
            s0 = s0+1;
            s1 = s0; 
        end
    end % for loop
    dr.ghxx = ghxx;
    dr.ghxu = ghxu;
    dr.ghuu = ghuu;
end

