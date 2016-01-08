function imcforecast(constrained_paths, constrained_vars, options_cond_fcst)
% Computes conditional forecasts.
%
% INPUTS
%  o constrained_paths    [double]      m*p array, where m is the number of constrained endogenous variables and p is the number of constrained periods. 
%  o constrained_vars     [char]        m*x array holding the names of the controlled endogenous variables. 
%  o options_cond_fcst    [structure]   containing the options. The fields are:
%                                                             + replic              [integer]   scalar, number of monte carlo simulations.
%                                                             + parameter_set       [char]      values of the estimated parameters:
%                                                                                               "posterior_mode", 
%                                                                                               "posterior_mean", 
%                                                                                               "posterior_median", 
%                                                                                               "prior_mode" or 
%                                                                                               "prior mean". 
%                                                                                   [double]     np*1 array, values of the estimated parameters.
%                                                             + controlled_varexo   [char]       m*x array, list of controlled exogenous variables.
%                                                             + conf_sig            [double]     scalar in [0,1], probability mass covered by the confidence bands. 
%
% OUTPUTS
%  None.
% 
% SPECIAL REQUIREMENTS
%  This routine has to be called after an estimation statement or an estimated_params block.
%
% REMARKS
%  [1] Results are stored in a structure which is saved in a mat file called conditional_forecasts.mat.
%  [2] Use the function plot_icforecast to plot the results. 

% Copyright (C) 2006-2010 Dynare Team
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

global options_ oo_ M_ bayestopt_

if isfield(options_cond_fcst,'parameter_set') || isempty(options_cond_fcst.parameter_set)
    options_cond_fcst.parameter_set = 'posterior_mode';
end

if ischar(options_cond_fcst.parameter_set)
    switch options_cond_fcst.parameter_set
      case 'posterior_mode'
        xparam = get_posterior_parameters('mode');
      case 'posterior_mean'
        xparam = get_posterior_parameters('mean');
      case 'posterior_median'
        xparam = get_posterior_parameters('median');
      case 'prior_mode'
        xparam = bayestopt_.p5(:);
      case 'prior_mean'
        xparam = bayestopt_.p1;
      otherwise
        disp('eval_likelihood:: If the input argument is a string, then it has to be equal to:')
        disp('                   ''posterior_mode'', ')
        disp('                   ''posterior_mean'', ')
        disp('                   ''posterior_median'', ')
        disp('                   ''prior_mode'' or')
        disp('                   ''prior_mean''.')
        error('imcforecast:: Wrong argument type!')
    end
else
    xparam = options_cond_fcst.parameter_set;
    if length(xparam)~=length(M_.params)
        error('imcforecast:: The dimension of the vector of parameters doesn''t match the number of estimated parameters!')
    end
end

if ~isfield(options_cond_fcst,'replic') || isempty(options_cond_fcst.replic)
    options_cond_fcst.replic = 5000;
end

if ~isfield(options_cond_fcst,'periods') || isempty(options_cond_fcst.periods)
    options_cond_fcst.periods = 40;
end

if ~isfield(options_cond_fcst,'conf_sig') || isempty(options_cond_fcst.conf_sig)
    options_cond_fcst.conf_sig = .8;
end

set_parameters(xparam);

n_varobs = size(options_.varobs,1);
rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);
options_ = set_default_option(options_,'nobs',size(rawdata,1)-options_.first_obs+1);
gend = options_.nobs;
rawdata = rawdata(options_.first_obs:options_.first_obs+gend-1,:);
% Transform the data.
if options_.loglinear
    if ~options_.logdata
        rawdata = log(rawdata);  
    end
end
% Test if the data set is real.
if ~isreal(rawdata)
    error('There are complex values in the data! Probably  a wrong transformation')
end
% Detrend the data.
options_.missing_data = any(any(isnan(rawdata)));
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
data = transpose(rawdata);
% Handle the missing observations.
[data_index,number_of_observations,no_more_missing_observations] = describe_missing_data(data,gend,n_varobs);
missing_value = ~(number_of_observations == gend*n_varobs);

[atT,innov,measurement_error,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(xparam,gend,data,data_index,number_of_observations);

trend = repmat(ys,1,options_cond_fcst.periods+1);
for i=1:M_.endo_nbr
    j = strmatch(deblank(M_.endo_names(i,:)),options_.varobs,'exact');
    if ~isempty(j)
        trend(i,:) = trend(i,:)+trend_coeff(j)*(gend+(0:options_cond_fcst.periods));
    end
end
trend = trend(oo_.dr.order_var,:);

InitState(:,1) = atT(:,end);
[T,R,ys,info] = dynare_resolve;

sQ = sqrt(M_.Sigma_e);

NumberOfStates = length(InitState);
FORCS1 = zeros(NumberOfStates,options_cond_fcst.periods+1,options_cond_fcst.replic);

for b=1:options_cond_fcst.replic
    FORCS1(:,1,b) = InitState;
end

EndoSize = M_.endo_nbr;
ExoSize = M_.exo_nbr;

n1 = size(constrained_vars,1);
n2 = size(options_cond_fcst.controlled_varexo,1);

if n1 ~= n2
    error(['imcforecast:: The number of constrained variables doesn''t match the number of controlled shocks'])
end

idx = [];
jdx = [];

for i = 1:n1
    idx = [idx ; oo_.dr.inv_order_var(strmatch(deblank(constrained_vars(i,:)),M_.endo_names,'exact'))];
    jdx = [jdx ; strmatch(deblank(options_cond_fcst.controlled_varexo(i,:)),M_.exo_names,'exact')];
end
mv = zeros(n1,NumberOfStates);
mu = zeros(ExoSize,n2);
for i=1:n1
    mv(i,idx(i)) = 1;
    mu(jdx(i),i) = 1;
end

if (size(constrained_paths,2) == 1);
    constrained_paths = constrained_paths*ones(1,cL);
else
    cL = size(constrained_paths,2);
end

constrained_paths = bsxfun(@minus,constrained_paths,trend(idx,1:cL));

randn('state',0);

for b=1:options_cond_fcst.replic
    shocks = sQ*randn(ExoSize,options_cond_fcst.periods);
    shocks(jdx,:) = zeros(length(jdx),options_cond_fcst.periods);
    FORCS1(:,:,b) = mcforecast3(cL,options_cond_fcst.periods,constrained_paths,shocks,FORCS1(:,:,b),T,R,mv, mu)+trend;
end

mFORCS1 = mean(FORCS1,3);

tt = (1-options_cond_fcst.conf_sig)/2;
t1 = round(options_cond_fcst.replic*tt);
t2 = round(options_cond_fcst.replic*(1-tt));

forecasts.controled_variables = constrained_vars;
forecasts.instruments = options_cond_fcst.controlled_varexo;

for i = 1:EndoSize
    eval(['forecasts.cond.mean.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = mFORCS1(i,:)'';']);
    tmp = sort(squeeze(FORCS1(i,:,:))');
    eval(['forecasts.cond.ci.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ...
          ' = [tmp(t1,:)'' ,tmp(t2,:)'' ]'';']);
end

clear FORCS1;

FORCS2 = zeros(NumberOfStates,options_cond_fcst.periods+1,options_cond_fcst.replic);
for b=1:options_cond_fcst.replic
    FORCS2(:,1,b) = InitState;
end

randn('state',0);

for b=1:options_cond_fcst.replic
    shocks = sQ*randn(ExoSize,options_cond_fcst.periods);
    shocks(jdx,:) = zeros(length(jdx),options_cond_fcst.periods);
    FORCS2(:,:,b) = mcforecast3(0,options_cond_fcst.periods,constrained_paths,shocks,FORCS2(:,:,b),T,R,mv, mu)+trend;
end

mFORCS2 = mean(FORCS2,3);

for i = 1:EndoSize
    eval(['forecasts.uncond.mean.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = mFORCS2(i,:)'';']);
    tmp = sort(squeeze(FORCS2(i,:,:))');
    eval(['forecasts.uncond.ci.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ...
          ' = [tmp(t1,:)'' ,tmp(t2,:)'' ]'';']);
end

save('conditional_forecasts.mat','forecasts');