function model_diagnostics(M_,options_,oo_)
% function model_diagnostics(M_,options_,oo_)
%   computes various diagnostics on the model 
% INPUTS
%   M_         [matlab structure] Definition of the model.           
%   options_   [matlab structure] Global options.
%   oo_        [matlab structure] Results 
%    
% OUTPUTS
%   none
%    
% ALGORITHM
%   ...
%    
% SPECIAL REQUIREMENTS
%   none.
%  

% Copyright (C) 1996-2011 Dynare Team
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

global jacob

endo_nbr = M_.endo_nbr;
endo_names = M_.endo_names;
lead_lag_incidence = M_.lead_lag_incidence;
maximum_lag = M_.maximum_lag;
maximum_lead = M_.maximum_lead;

%
% missing variables at the current period
%
k = find(lead_lag_incidence(maximum_lag+1,:)==0);
if ~isempty(k)
    disp(['The following endogenous variables aren''t present at ' ...
          'the current period in the model:'])
    for i=1:length(k)
        disp(endo_names(k(i),:))
    end
end

%
% check steady state
%
info = 0;

it_ = M_.maximum_lag + 1 ;

if M_.exo_nbr == 0
    oo_.exo_steady_state = [] ;
end

% check if ys is steady state
tempex = oo_.exo_simul;
oo_.exo_simul = repmat(oo_.exo_steady_state',M_.maximum_lag+M_.maximum_lead+1,1);
if M_.exo_det_nbr > 0 
    tempexdet = oo_.exo_det_simul;
    oo_.exo_det_simul = repmat(oo_.exo_det_steady_state',M_.maximum_lag+M_.maximum_lead+1,1);
end
dr.ys = oo_.steady_state;
check1 = 0;
% testing for steadystate file
fh = str2func([M_.fname '_static']);
if options_.steadystate_flag
    [ys,check1] = feval([M_.fname '_steadystate'],dr.ys,...
                        [oo_.exo_steady_state; oo_.exo_det_steady_state]);
    M_.params = evalin('base','M_.params;');
    if size(ys,1) < M_.endo_nbr 
        if length(M_.aux_vars) > 0
            ys = add_auxiliary_variables_to_steadystate(ys,M_.aux_vars,...
                                                        M_.fname,...
                                                        oo_.exo_steady_state,...
                                                        oo_.exo_det_steady_state,...
                                                        options_.bytecode);
        else
            error([M_.fname '_steadystate.m doesn''t match the model']);
        end
    end
    dr.ys = ys;
else
    % testing if ys isn't a steady state or if we aren't computing Ramsey policy
    if  options_.ramsey_policy == 0
        if options_.linear == 0
            % nonlinear models
            if max(abs(feval(fh,dr.ys,[oo_.exo_steady_state; ...
                                    oo_.exo_det_steady_state], M_.params))) > options_.dynatol
                [ys,check1] = dynare_solve(fh,dr.ys,1,...
                                           [oo_.exo_steady_state; ...
                                    oo_.exo_det_steady_state], ...
                                           M_.params);
                if ~check1
                    dr.ys = ys;
                end
            end
        else
            % linear models
            [fvec,jacob] = feval(fh,dr.ys,[oo_.exo_steady_state;...
                                oo_.exo_det_steady_state], M_.params);
            if max(abs(fvec)) > 1e-12
                dr.ys = dr.ys-jacob\fvec;
            end
        end
    end
end
% testing for problem
if check1
    disp('model diagnostic can''t obtain the steady state')
end

if ~isreal(dr.ys)
    disp(['model diagnostic obtains a steady state with complex ' ...
          'numbers'])
    return
end

%
% singular Jacobian of static model
%

[res,jacob]=feval(fh,dr.ys,[oo_.exo_steady_state; oo_.exo_det_steady_state], ...
                  M_.params);
rank_jacob = rank(jacob);
if rank_jacob < endo_nbr
    disp(['model_diagnostic: the Jacobian of the static model is ' ...
          'singular'])
    disp(['there is ' num2str(endo_nbr-rank_jacob) ...
          ' colinear relationships between the variables and the equations'])
    ncol = null(jacob);
    n_rel = size(ncol,2);
    for i = 1:n_rel
        if n_rel  > 1
            disp(['Relation ' int2str(i)])
        end
        disp('Colinear variables:')
        for j=1:10
            k = find(abs(ncol(:,i)) > 10^-j);
            if max(abs(jacob(:,k)*ncol(k,i))) < 1e-6
                break
            end
        end
        disp(endo_names(k,:))
    end
    neq = null(jacob');
    n_rel = size(neq,2);
    for i = 1:n_rel
        if n_rel  > 1
            disp(['Relation ' int2str(i)])
        end
        disp('Colinear equations')
        for j=1:10
            k = find(abs(neq(:,i)) > 10^-j);
            if max(abs(jacob(k,:)'*neq(k,i))) < 1e-6
                break
            end
        end
        disp(k')
    end
end

