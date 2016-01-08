function osr1(i_params,i_var,weights)

% Copyright (C) 2005-2011 Dynare Team
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

global M_ oo_ options_ it_

klen = M_.maximum_lag + M_.maximum_lead + 1;
iyv = M_.lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = M_.maximum_lag + 1 ;


if M_.exo_nbr == 0
    oo_.exo_steady_state = [] ;
end

if ~ M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0
    error ('OSR: Error in model specification: some variables don''t appear as current') ;
end

if M_.maximum_lead == 0
    error ('Backward or static model: no point in using OSR') ;
end

exe =zeros(M_.exo_nbr,1);

dr = set_state_space(oo_.dr,M_);

% check if ys is steady state
if exist([M_.fname '_steadystate'])
    [ys,check1] = feval([M_.fname '_steadystate'],oo_.steady_state,...
                        [oo_.exo_steady_state; oo_.exo_det_steady_state]);
    if size(ys,1) < M_.endo_nbr 
        if length(M_.aux_vars) > 0
            ys = add_auxiliary_variables_to_steadystate(ys,M_.aux_vars,...
                                                        M_.fname,...
                                                        oo_.exo_steady_state,...
                                                        oo_.exo_det_steady_state,M_.params,...
                                                        options_.bytecode);
        else
            error([M_.fname '_steadystate.m doesn''t match the model']);
        end
    end
    dr.ys = ys;
else
    % testing if ys isn't a steady state or if we aren't computing Ramsey policy
    fh = str2func([M_.fname '_static']);
    if max(abs(feval(fh,oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state], M_.params))) ...
            > options_.dynatol && options_.ramsey_policy == 0
        if options_.linear == 0
            % nonlinear models
            [dr.ys,check1] = dynare_solve(fh,dr.ys,options_.jacobian_flag,...
                                          [oo_.exo_steady_state; ...
                                oo_.exo_det_steady_state], M_.params);
        else
            % linear models
            [fvec,jacob] = feval(fh,oo_.steady_state,[oo_.exo_steady_state;...
                                oo_.exo_det_steady_state], M_.params);
            dr.ys = oo_.steady_state-jacob\fvec;
        end
    end
end
oo_.dr = dr;
% $$$   if max(abs(feval(fh, oo_.steady_state, exe, M_.params))) > options_.dynatol
% $$$     [oo_.dr.ys, check] = dynare_solve([M_.fname '_static'], oo_.steady_state, 1, exe, M_.params);
% $$$     if check
% $$$       error('OLR: convergence problem in DYNARE_SOLVE')
% $$$     end
% $$$   else
% $$$     oo_.dr.ys = oo_.steady_state;
% $$$   end


np = size(i_params,1);
t0 = M_.params(i_params);
inv_order_var = oo_.dr.inv_order_var;

H0 = 1e-4*eye(np);
crit = 1e-7;
nit = 1000;
verbose = 2;

[f,p]=csminwel1('osr_obj',t0,H0,[],crit,nit,options_.gradient_method,options_.gradient_epsilon,i_params,...
                inv_order_var(i_var),weights(i_var,i_var));

%  options = optimset('fminunc');
%  options = optimset('display','iter');
%  [p,f]=fminunc(@osr_obj,t0,options,i_params,...
%               inv_order_var(i_var),weights(i_var,i_var));



disp('')
disp('OPTIMAL VALUE OF THE PARAMETERS:')
disp('')
for i=1:np
    disp(sprintf('%16s %16.6g\n',M_.param_names(i_params(i),:),p(i)))
end
disp(sprintf('Objective function : %16.6g\n',f));
disp(' ')
oo_.dr=resol(oo_.steady_state,0);

% 05/10/03 MJ modified to work with osr.m and give full report