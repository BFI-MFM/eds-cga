function [data,rawdata,xparam1,data_info]=dynare_estimation_init(var_list_, dname, gsa_flag)

% function dynare_estimation_init(var_list_, gsa_flag)
% preforms initialization tasks before estimation or
% global sensitivity analysis
%  
% INPUTS
%   var_list_:  selected endogenous variables vector
%   dname:      alternative directory name
%   gsa_flag:   flag for GSA operation (optional)
%  
% OUTPUTS
%   data:    data after required transformation
%   rawdata:  data as in the data file
%   xparam1:    initial value of estimated parameters as returned by
%               set_prior()
%
% SPECIAL REQUIREMENTS
%   none

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

global M_ options_ oo_ estim_params_ bayestopt_

if nargin < 3 || isempty(gsa_flag)
    gsa_flag = 0;
else
    %% Decide if a DSGE or DSGE-VAR has to be estimated.
    if ~isempty(strmatch('dsge_prior_weight',M_.param_names))
        options_.dsge_var = 1;
    end
    var_list_ = check_list_of_variables(options_, M_, var_list_);
    options_.varlist = var_list_;
end

options_.lgyidx2varobs = zeros(size(M_.endo_names,1),1);
for i = 1:size(M_.endo_names,1)
    tmp = strmatch(deblank(M_.endo_names(i,:)),options_.varobs,'exact');
    if ~isempty(tmp)
        if length(tmp)>1
            disp(' ')
            error(['Multiple declarations of ' deblank(M_.endo_names(i,:)) ' as an observed variable is not allowed!'])
        end
        options_.lgyidx2varobs(i,1) = tmp;
    end
end

%% Set the order of approximation to one (if needed).
if options_.order > 1
    if ~exist('particle','dir')
        disp('This version of Dynare cannot estimate non linearized models!')
        disp('Set "order" equal to 1.')
        disp(' ')
        options_.order = 1;
    end
end

%% Set options_.lik_init equal to 3 if diffuse filter is used.
if (options_.diffuse_filter==1) && (options_.lik_init==1)
    options_.lik_init = 3;
end

%% If options_.lik_init == 1
%%  set by default options_.qz_criterium to 1-1e-6 
%%  and check options_.qz_criterium < 1-eps if options_.lik_init == 1
%% Else set by default options_.qz_criterium to 1+1e-6
if options_.lik_init == 1
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1-1e-6;
    elseif options_.qz_criterium > 1-eps
        error(['estimation: option qz_criterium is too large for estimating ' ...
               'a stationary model. If your model contains unit roots, use ' ...
               'option diffuse_filter'])
    end
elseif isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

%% If the data are prefiltered then there must not be constants in the
%% measurement equation of the DSGE model or in the DSGE-VAR model.
if options_.prefilter == 1
    options_.noconstant = 1;
end

%% Set options related to filtered variables.
if options_.filtered_vars ~= 0 && isempty(options_.filter_step_ahead), 
    options_.filter_step_ahead = 1;
end

if options_.filter_step_ahead ~= 0
    options_.nk = max(options_.filter_step_ahead);
end

%% Set the name of the directory where (intermediary) results will be saved.
if nargin>1
    M_.dname = dname;
else
    M_.dname = M_.fname; 
end

%% Set the number of observed variables.
n_varobs = size(options_.varobs,1);

%% Set priors over the estimated parameters.
if ~isempty(estim_params_)
    [xparam1,estim_params_,bayestopt_,lb,ub,M_] = set_prior(estim_params_,M_,options_);
    if any(bayestopt_.pshape > 0)
        % Plot prior densities.
        if options_.plot_priors
            plot_priors(bayestopt_,M_,options_)
        end
        % Set prior bounds
        bounds = prior_bounds(bayestopt_);
        bounds(:,1)=max(bounds(:,1),lb);
        bounds(:,2)=min(bounds(:,2),ub);
    else
        % No priors are declared so Dynare will estimate the model by
        % maximum likelihood with inequality constraints for the parameters.
        options_.mh_replic = 0;% No metropolis.
        bounds(:,1) = lb;
        bounds(:,2) = ub;
    end
    % Test if initial values of the estimated parameters are all between
    % the prior lower and upper bounds.
    if any(xparam1 < bounds(:,1)) || any(xparam1 > bounds(:,2))
        find(xparam1 < bounds(:,1))
        find(xparam1 > bounds(:,2))
        error('Initial parameter values are outside parameter bounds')
    end
    lb = bounds(:,1);
    ub = bounds(:,2);
    bayestopt_.lb = lb;
    bayestopt_.ub = ub;
else% If estim_params_ is empty...
    xparam1 = [];
    bayestopt_.lb = [];
    bayestopt_.ub = [];
    bayestopt_.jscale = [];
    bayestopt_.pshape = [];
    bayestopt_.p1 = [];
    bayestopt_.p2 = [];
    bayestopt_.p3 = [];
    bayestopt_.p4 = [];
    bayestopt_.p5 = [];
    bayestopt_.p6 = [];
    bayestopt_.p7 = [];
    estim_params_.nvx = 0;
    estim_params_.nvn = 0;
    estim_params_.ncx = 0;
    estim_params_.ncn = 0;
    estim_params_.np = 0;
end

%% Is there a linear trend in the measurement equation?
if ~isfield(options_,'trend_coeffs') % No!
    bayestopt_.with_trend = 0;
else% Yes!
    bayestopt_.with_trend = 1;
    bayestopt_.trend_coeff = {};
    trend_coeffs = options_.trend_coeffs;
    nt = length(trend_coeffs);
    for i=1:n_varobs
        if i > length(trend_coeffs)
            bayestopt_.trend_coeff{i} = '0';
        else
            bayestopt_.trend_coeff{i} = trend_coeffs{i};
        end
    end
end

%% Set the "size" of penalty.
bayestopt_.penalty = 1e8; 

%% Get informations about the variables of the model.
dr = set_state_space(oo_.dr,M_);
oo_.dr = dr;
nstatic = dr.nstatic;          % Number of static variables. 
npred = dr.npred;              % Number of predetermined variables.
nspred = dr.nspred;            % Number of predetermined variables in the state equation.

%% Test if observed variables are declared.
if isempty(options_.varobs)
    error('VAROBS is missing')
end

%% Setting resticted state space (observed + predetermined variables)
var_obs_index = [];
k1 = [];
for i=1:n_varobs
    var_obs_index = [var_obs_index strmatch(deblank(options_.varobs(i,:)),M_.endo_names(dr.order_var,:),'exact')];
    k1 = [k1 strmatch(deblank(options_.varobs(i,:)),M_.endo_names, 'exact')];
end
% Define union of observed and state variables
k2 = union(var_obs_index',[dr.nstatic+1:dr.nstatic+dr.npred]', 'rows');
% Set restrict_state to postion of observed + state variables in expanded state vector.
oo_.dr.restrict_var_list = k2;
% set mf0 to positions of state variables in restricted state vector for likelihood computation.
[junk,bayestopt_.mf0] = ismember([dr.nstatic+1:dr.nstatic+dr.npred]',k2);
% Set mf1 to positions of observed variables in restricted state vector for likelihood computation.
[junk,bayestopt_.mf1] = ismember(var_obs_index,k2); 
% Set mf2 to positions of observed variables in expanded state vector for filtering and smoothing.
bayestopt_.mf2  = var_obs_index;
bayestopt_.mfys = k1;

[junk,ic] = intersect(k2,nstatic+(1:npred)');
oo_.dr.restrict_columns = [ic; length(k2)+(1:nspred-npred)'];

k3 = [];
if options_.selected_variables_only
    for i=1:size(var_list_,1)
        k3 = [k3; strmatch(var_list_(i,:),M_.endo_names(dr.order_var,:), ...
                           'exact')];
    end
else
    k3 = (1:M_.endo_nbr)';
end
bayestopt_.smoother_var_list = union(k2,k3);
[junk,bayestopt_.smoother_saved_var_list] = intersect(k3,bayestopt_.smoother_var_list(:));
[junk,ic] = intersect(bayestopt_.smoother_var_list,nstatic+(1:npred)');
bayestopt_.smoother_restrict_columns = ic;
[junk,bayestopt_.smoother_mf] = ismember(var_obs_index, ...
                                         bayestopt_.smoother_var_list);

%% Initialization with unit-root variables.
if ~isempty(options_.unit_root_vars)
    n_ur = size(options_.unit_root_vars,1);
    i_ur = zeros(n_ur,1);
    for i=1:n_ur
        i1 = strmatch(deblank(options_.unit_root_vars(i,:)),M_.endo_names(dr.order_var,:),'exact');
        if isempty(i1)
            error('Undeclared variable in unit_root_vars statement')
        end
        i_ur(i) = i1;
    end
    bayestopt_.var_list_stationary = setdiff((1:M_.endo_nbr)',i_ur);
    [junk,bayestopt_.restrict_var_list_nonstationary] = ...
        intersect(oo_.dr.restrict_var_list,i_ur);
    bayestopt_.restrict_var_list_stationary = ...
        setdiff((1:length(oo_.dr.restrict_var_list))', ...
                bayestopt_.restrict_var_list_nonstationary);
    if M_.maximum_lag > 1
        l1 = flipud([cumsum(M_.lead_lag_incidence(1:M_.maximum_lag-1,dr.order_var),1);ones(1,M_.endo_nbr)]);
        l2 = l1(:,oo_.dr.restrict_var_list);
        il2 = find(l2' > 0);
        l2(il2) = (1:length(il2))';
        bayestopt_.restrict_var_list_stationary = ...
            nonzeros(l2(:,bayestopt_.restrict_var_list_stationary)); 
        bayestopt_.restrict_var_list_nonstationary = ...
            nonzeros(l2(:,bayestopt_.restrict_var_list_nonstationary)); 
    end
    options_.lik_init = 3;
end % if ~isempty(options_.unit_root_vars)

%% Test if the data file is declared.
if isempty(options_.datafile)
    if gsa_flag
        data = [];
        rawdata = [];
        data_info = [];
        return
    else
        error('datafile option is missing')
    end     
end

%% If jscale isn't specified for an estimated parameter, use global option options_.jscale, set to 0.2, by default.
k = find(isnan(bayestopt_.jscale));
bayestopt_.jscale(k) = options_.mh_jscale;

%% Load and transform data.
rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);
% Set the number of observations (nobs) and build a subsample between first_obs and nobs.
options_ = set_default_option(options_,'nobs',size(rawdata,1)-options_.first_obs+1);
gend = options_.nobs;
rawdata = rawdata(options_.first_obs:options_.first_obs+gend-1,:);
% Take the log of the variables if needed
if options_.loglinear      % If the model is log-linearized...
    if ~options_.logdata   % and if the data are not in logs, then...
        rawdata = log(rawdata);  
    end
end
% Test if the observations are real numbers. 
if ~isreal(rawdata)
    error('There are complex values in the data! Probably  a wrong transformation')
end
% Test for missing observations.
options_.missing_data = any(any(isnan(rawdata)));
% Prefilter the data if needed.
if options_.prefilter == 1
    if options_.missing_data
        bayestopt_.mean_varobs = zeros(n_varobs,1);
        for variable=1:n_varobs
            rdx = find(~isnan(rawdata(:,variable)));
            m = mean(rawdata(rdx,variable));
            rawdata(rdx,variable) = rawdata(rdx,variable)-m;
            bayestopt_.mean_varobs(variable) = m;
        end
    else
        bayestopt_.mean_varobs = mean(rawdata,1)';
        rawdata = rawdata-repmat(bayestopt_.mean_varobs',gend,1);
    end
end
% Transpose the dataset array.
data = transpose(rawdata);

if nargout>3
    %% Compute the steady state: 
    if options_.steadystate_flag% if the *_steadystate.m file is provided.
        [ys,tchek] = feval([M_.fname '_steadystate'],...
                           [zeros(M_.exo_nbr,1);...
                            oo_.exo_det_steady_state]);
        if size(ys,1) < M_.endo_nbr 
            if length(M_.aux_vars) > 0
                ys = add_auxiliary_variables_to_steadystate(ys,M_.aux_vars,...
                                                            M_.fname,...
                                                            zeros(M_.exo_nbr,1),...
                                                            oo_.exo_det_steady_state,...
                                                            M_.params,...
                                                            options_.bytecode);
            else
                error([M_.fname '_steadystate.m doesn''t match the model']);
            end
        end
        oo_.steady_state = ys;
    else% if the steady state file is not provided.
        [dd,info] = resol(oo_.steady_state,0);
        oo_.steady_state = dd.ys; clear('dd');
    end
    if all(abs(oo_.steady_state(bayestopt_.mfys))<1e-9)
        options_.noconstant = 1;
    else
        options_.noconstant = 0;
    end
    
    [data_index,number_of_observations,no_more_missing_observations] = describe_missing_data(data,gend,n_varobs);
    missing_value = ~(number_of_observations == gend*n_varobs);
    
%     initial_estimation_checks(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations);

    data_info.gend = gend;
    data_info.data = data;
    data_info.data_index = data_index;
    data_info.number_of_observations = number_of_observations;
    data_info.no_more_missing_observations = no_more_missing_observations;
    data_info.missing_value = missing_value;
end
