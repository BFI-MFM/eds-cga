function innovation_paths = reversed_extended_path(controlled_variable_names, control_innovation_names, dataset)
% Inversion of the extended path simulation approach. This routine computes the innovations needed to
% reproduce the time path of a subset of endogenous variables. The initial condition is teh deterministic
% steady state.   
%    
% INPUTS    
%  o controlled_variable_names        [string]    n*1 matlab's cell. 
%  o control_innovation_names         [string]    n*1 matlab's cell.  
%  o dataset                          [structure]
% OUTPUTS
%  o innovations                      [double]  n*T matrix.
%    
% ALGORITHM
%  
% SPECIAL REQUIREMENTS

% Copyright (C) 2010 Dynare Team.
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

global M_ oo_ options_   

%% Initialization

% Load data.
eval(dataset.name);
dataset.data = [];
for v = 1:dataset.number_of_observed_variables
    eval(['dataset.data = [ dataset.data , ' dataset.variables(v,:) ' ];'])
end
data = dataset.data(dataset.first_observation:dataset.first_observation+dataset.number_of_observations,:);

% Compute the deterministic steady state.
steady_;

%  Compute the first order perturbation reduced form.
old_options_order = options_.order; options_.order = 1;
[oo_.dr,info]  = resol(oo_.steady_state,0);
options_.order = old_options_order;

% Set various options.
options_.periods = 100;

% Set-up oo_.exo_simul.
make_ex_; 

% Set-up oo_.endo_simul.
make_y_;

% Get indices of the controlled endogenous variables in endo_simul.
n  = length(controlled_variable_names);
iy = NaN(n,1);
for k=1:n
    iy(k) = strmatch(controlled_variable_names{k},M_.endo_names,'exact');
end

% Get indices of the controlled endogenous variables in dataset.
iy_ = NaN(n,1);
for k=1:n
    iy_(k) = strmatch(controlled_variable_names{k},dataset.variables,'exact');
end

% Get indices of the control innovations in exo_simul.
ix = NaN(n,1);
for k=1:n
    ix(k) = strmatch(control_innovation_names{k},M_.exo_names,'exact');
end

% Get the length of the sample.
T = size(data,1);

% Output initialization.
innovation_paths = zeros(n,T);

% Initialization of the perfect foresight model solver.
perfect_foresight_simulation();

% Set options for fsolve.
options = optimset('MaxIter',10000,'Display','Iter');

%% Call fsolve recursively
for t=1:T
    x0 = zeros(n,1);
    y_target  = transpose(data(t,iy_));
    total_variation = y_target-transpose(oo_.endo_simul(t+M_.maximum_lag,iy));
    for i=1:100
        [t,i]
        y = transpose(oo_.endo_simul(t+M_.maximum_lag,iy)) + (i/100)*y_target 
        [tmp,fval,exitflag] = fsolve('ep_residuals', x0, options, y, ix, iy, oo_.steady_state, oo_.dr, M_.maximum_lag, M_.endo_nbr);
    end
    if exitflag==1
        innovation_paths(:,t) = tmp;
    end
    % Update endo_simul.
    oo_.endo_simul(:,1:end-1) = oo_.endo_simul(:,2:end); 
    oo_.endo_simul(:,end) = oo_.steady_state;
end