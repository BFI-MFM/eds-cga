function global_initialization()
%function global_initialization()
% initializes global variables and options for DYNARE
%
% INPUTS
%    none
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

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

global oo_ M_ options_ estim_params_

estim_params_ = [];

options_.console_mode = 0;

options_.terminal_condition = 0;
options_.rplottype = 0;
options_.smpl = 0;
options_.dynatol = 0.00001;
options_.maxit_ = 10;
options_.slowc = 1;
options_.timing = 0;
options_.gstep = 1e-2;
options_.scalv = 1;
options_.debug = 0;
options_.initval_file = 0;
options_.Schur_vec_tol = 1e-11; % used to find nonstationary variables
                                % in Schur decomposition of the
                                % transition matrix
options_.qz_criterium = [];
options_.lyapunov_complex_threshold = 1e-15;
options_.solve_tolf = eps^(1/3);
options_.solve_tolx = eps^(2/3);
options_.solve_maxit = 500;
options_.deterministic_simulation_initialization = 0;

% Default number of threads for parallelized mex files.  
options_.threads.kronecker.A_times_B_kronecker_C = 1;
options_.threads.kronecker.sparse_hessian_times_B_kronecker_C = 1;

% steady state file
if exist([M_.fname '_steadystate.m'],'file')
    options_.steadystate_flag = 1;
else
    options_.steadystate_flag = 0;
end
options_.steadystate_partial = [];

% subset of the estimated deep parameters
options_.ParamSubSet = 'None';

% bvar-dsge
options_.dsge_var = 0;
options_.dsge_varlag = 4;

% Optimization algorithm [6] gmhmaxlik
options_.Opt6Iter = 2;
options_.Opt6Numb = 5000;

% Graphics
options_.graphics.nrows = 3;
options_.graphics.ncols = 3;
options_.graphics.line_types = {'b-'};
options_.graphics.line_width = 1;
options_.nograph = 0;
options_.XTick = [];
options_.XTickLabel = [];

% IRFs & other stoch_simul output
options_.irf = 40;
options_.relative_irf = 0;
options_.ar = 5;
options_.hp_filter = 0;
options_.hp_ngrid = 512;
options_.nomoments = 0;
options_.nocorr = 0;
options_.periods = 0;
options_.noprint = 0;
options_.SpectralDensity = 0;

% TeX output
options_.TeX = 0;

% Exel
options_.xls_sheet = '';
options_.xls_range = '';

% Prior draws
options_.forecast = 0;

% Model
options_.linear = 0;

% Deterministic simulation
options_.stack_solve_algo = 0;
options_.markowitz = 0.5;
options_.minimal_solving_periods = 1;

% Solution
options_.order = 2;
options_.pruning = 0;
options_.solve_algo = 2;
options_.linear = 0;
options_.replic = 50;
options_.drop = 100;
% if mjdgges.dll (or .mexw32 or ....) doesn't exist, matlab/qz is added to the path. 
% There exists now qz/mjdgges.m that contains the calls to the old Sims code 
% Hence, if mjdgges.m is visible exist(...)==2, 
% this means that the DLL isn't avaiable and use_qzdiv is set to 1
if exist('mjdgges')==2
    options_.use_qzdiv = 1;
else
    options_.use_qzdiv = 0;
end
options_.aim_solver = 0; % i.e. by default do not use G.Anderson's AIM solver, use mjdgges instead
options_.k_order_solver=0; % by default do not use k_order_perturbation but mjdgges
options_.partial_information = 0;
options_.ACES_solver = 0;
options_.conditional_variance_decomposition = [];

% Ramsey policy
options_.planner_discount = 1.0;
options_.ramsey_policy = 0;
options_.timeless = 0;

% estimation
options_.Harvey_scale_factor = 10;
options_.MaxNumberOfBytes = 1e6;
options_.MaximumNumberOfMegaBytes = 111;
options_.PosteriorSampleSize = 1000;
options_.bayesian_irf = 0;
options_.bayesian_th_moments = 0;
options_.diffuse_filter = 0;
options_.filter_step_ahead = [];
options_.filtered_vars = 0;
options_.first_obs = 1;
options_.kalman_algo = 0;
options_.kalman_tol = 1e-10;
options_.riccati_tol = 1e-6;
options_.lik_algo = 1;
options_.lik_init = 1;
options_.load_mh_file = 0;
options_.logdata = 0;
options_.loglinear = 0;
options_.mh_conf_sig = 0.90;
options_.prior_interval = 0.90;
options_.mh_drop = 0.5;
options_.mh_jscale = 0.2;
options_.mh_init_scale = 2*options_.mh_jscale;
options_.mh_mode = 1;
options_.mh_nblck = 2;
options_.mh_recover = 0;
options_.mh_replic = 20000;
options_.mode_check = 0;
options_.mode_check_nolik = 0;
options_.mode_compute = 4;
options_.mode_file = '';
options_.moments_varendo = 0;
options_.nk = 1;
options_.noconstant = 0;
options_.nodiagnostic = 0;
options_.mh_posterior_mode_estimation = 0;
options_.prefilter = 0;
options_.presample = 0;
options_.prior_trunc = 1e-10;
options_.smoother = 0;
options_.student_degrees_of_freedom = 3;
options_.subdraws = [];
options_.unit_root_vars = [];
options_.use_mh_covariance_matrix = 0;
options_.gradient_method = 2;
options_.gradient_epsilon = 1e-6;
options_.posterior_sampling_method = 'random_walk_metropolis_hastings';
options_.proposal_distribution = 'rand_multivariate_normal';
options_.student_degrees_of_freedom = 3;
options_.trace_plot_ma = 200;
options_.mh_autocorrelation_function_size = 30;
options_.plot_priors = 1;
options_.cova_compute = 1;
options_.parallel = 0;
options_.parallel_info.leaveSlaveOpen = 0;
options_.parallel_info.RemoteTmpFolder = '';
options_.number_of_grid_points_for_kde = 2^9;
quarter = 1;
years = [1 2 3 4 5 10 20 30 40 50];
options_.conditional_variance_decomposition_dates = zeros(1,length(years));
for i=1:length(years)
    options_.conditional_variance_decomposition_dates(i) = ...
        (years(i)-1)*4+quarter;
end
options_.filter_covariance = 0;
options_.filter_decomposition = 0;
options_.selected_variables_only = 0;
% Misc
options_.conf_sig = 0.6;
oo_.exo_simul = [];
oo_.endo_simul = [];
oo_.dr = [];
oo_.exo_steady_state = [];
oo_.exo_det_steady_state = [];
oo_.exo_det_simul = [];

M_.params = [];

% BVAR
M_.bvar = [];

% homotopy
options_.homotopy_mode = 0;
options_.homotopy_steps = 1;

% prior analysis
options_.prior_mc = 20000;
options_.prior_analysis_endo_var_list = [];

% did model undergo block decomposition + minimum feedback set computation ?
options_.block = 0;

% model evaluated using a compiled MEX
options_.use_dll = 0;

% model evaluated using bytecode.dll
options_.bytecode = 0;

% dates for historical time series
options_.initial_date.freq = 1;
options_.initial_date.period = 1;
options_.initial_date.subperiod = 0;

% SWZ SBVAR
options_.ms.freq = 1;
options_.ms.initial_subperiod = 1;
options_.ms.final_subperiod=4;
options_.ms.log_var = [ ];
options_.ms.forecast = 1;
options_.ms.nlags = 1;
options_.ms.cross_restrictions = 0;
options_.ms.contemp_reduced_form = 0;
options_.ms.real_pseudo_forecast = 0;
options_.ms.bayesian_prior = 1;
options_.ms.dummy_obs = 0;
options_.ms.ncsk = 0;
options_.ms.nstd = 6;
options_.ms.ninv = 1000;
options_.ms.indxparr = 1;
options_.ms.indxovr = 0;
options_.ms.aband = 1;
options_.ms.indxap = 1;
options_.ms.apband = 1;
options_.ms.indximf = 1;
options_.ms.imfband = 1;
options_.ms.indxfore = 0;
options_.ms.foreband = 0;
options_.ms.indxgforhat = 1;
options_.ms.indxgimfhat = 1;
options_.ms.indxestima = 1;
options_.ms.indxgdls = 1;
options_.ms.cms =0;
options_.ms.ncms = 0;
options_.ms.eq_cms = 1;
options_.ms.cnum = 0;
options_.ms.banact = 1;
options_.ms.nstates = 2;
options_.ms.indxscalesstates = 0;
options_.ms.alpha = 1.0;
options_.ms.beta = 1.0;
options_.ms.gsig2_lmd = 50^2;
options_.ms.gsig2_lmdm = 50^2;
options_.ms.q_diag = 0.85;
options_.ms.flat_prior = 0;
options_.ms.create_initialization_file = 1;
options_.ms.estimate_msmodel = 1;
options_.ms.compute_mdd = 1;
options_.ms.compute_probabilities = 1;
options_.ms.print_draws = 1;
options_.ms.n_draws=1000;
options_.ms.thinning_factor=1;
options_.ms.proposal_draws = 100000;
options_.ms.lower_cholesky = 0;
options_.ms.upper_cholesky = 0;
options_.ms.Qi = [];
options_.ms.Ri = [];
options_.ms.draws_nbr_burn_in_1 = 30000;
options_.ms.draws_nbr_burn_in_2 = 10000;
options_.ms.draws_nbr_mean_var_estimate = 200000;
options_.ms.draws_nbr_modified_harmonic_mean = 1000000;
options_.ms.thinning_factor = 1;
options_.ms.dirichlet_scale = [1.0 1.5 2.0];

% Shock decomposition
options_.parameter_set = [];

% initialize persistent variables in priordens()
priordens([],[],[],[],[],[],1);

% Set dynare random generator and seed.
set_dynare_seed('default');

