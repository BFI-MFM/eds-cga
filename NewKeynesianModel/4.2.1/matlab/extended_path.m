function time_series = extended_path(initial_conditions,sample_size,init)
% Stochastic simulation of a non linear DSGE model using the Extended Path method (Fair and Taylor 1983). A time
% series of size T  is obtained by solving T perfect foresight models. 
%    
% INPUTS
%  o initial_conditions     [double]    m*nlags array, where m is the number of endogenous variables in the model and
%                                       nlags is the maximum number of lags.
%  o sample_size            [integer]   scalar, size of the sample to be simulated.
%  o init                   [integer]   scalar, method of initialization of the perfect foresight equilibrium paths
%                                                    init=0  previous solution is used,
%                                                    init=1  a path generated with the first order reduced form is used.
%                                                    init=2  mix of cases 0 and 1. 
%   
% OUTPUTS
%  o time_series            [double]    m*sample_size array, the simulations.
%    
% ALGORITHM
%  
% SPECIAL REQUIREMENTS

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
global M_ oo_ options_ 

% Set default initial conditions.
if isempty(initial_conditions) 
    initial_conditions = repmat(oo_.steady_state,1,M_.maximum_lag); 
end

% Set default value for the last input argument
if nargin<3
    init = 0;
end

% Set the number of periods for the deterministic solver.
%options_.periods = 40;

% Initialize the exogenous variables.
make_ex_; 

% Initialize the endogenous variables.
make_y_;

% Compute the first order reduced form if needed.
if init
    oldopt = options_;
    options_.order = 1;
    [dr,info]=resol(oo_.steady_state,0);
    oo_.dr = dr;
    options_ = oldopt;
    if init==2
        lambda = .8;
    end
end

% Initialize the output array.
time_series = NaN(M_.endo_nbr,sample_size+1); 

% Set the covariance matrix of the structural innovations.
variances = diag(M_.Sigma_e); 
positive_var_indx = find(variances>0); 
covariance_matrix = M_.Sigma_e(positive_var_indx,positive_var_indx); 
number_of_structural_innovations = length(covariance_matrix); 
covariance_matrix_upper_cholesky = chol(covariance_matrix); 

tdx = M_.maximum_lag+1; 
norme = 0;

% Set verbose option
verbose = 0;

t = 0;
new_draw = 1;

perfect_foresight_simulation();

while (t<=sample_size)
    t = t+1;
    if new_draw
        gaussian_draw = randn(1,number_of_structural_innovations);
    else
        gaussian_draw = .5*gaussian_draw ;
        new_draw = 1;
    end
    shocks = exp(gaussian_draw*covariance_matrix_upper_cholesky-.5*variances(positive_var_indx)');
    oo_.exo_simul(tdx,positive_var_indx) = shocks;
    if init
        % Compute first order solution.
        exogenous_variables = zeros(size(oo_.exo_simul));
        exogenous_variables(tdx,positive_var_indx) = log(shocks);
        initial_path = simult_(oo_.steady_state,dr,exogenous_variables,1);
        if init==1
            oo_.endo_simul = initial_path(:,1:end-1);
        else
            oo_.endo_simul = initial_path(:,1:end-1)*lambda + oo_.endo_simul*(1-lambda);  
        end
    end
    if init
        info = perfect_foresight_simulation(dr,oo_.steady_state);
    else
        info = perfect_foresight_simulation;
    end
    time = info.time;
    if verbose
        [t,options_.periods]
        info
        info.iterations
    end
    if ~info.convergence
        INFO = homotopic_steps(tdx,positive_var_indx,shocks,norme,.5,init,0);
        if verbose
            norme
            INFO
        end
        if ~isstruct(INFO) && isnan(INFO)
            t = t-1;
            new_draw = 0;
        else
            info = INFO;
        end
    else
        norme = sqrt(sum((shocks-1).^2,2));
    end
    %if ~info.convergence
    %    error('I am not able to simulate this model!')
    %end
    if new_draw
        info.time = info.time+time;
        time_series(:,t+1) = oo_.endo_simul(:,tdx);
        oo_.endo_simul(:,1:end-1) = oo_.endo_simul(:,2:end); 
        oo_.endo_simul(:,end) = oo_.steady_state;
    end
end