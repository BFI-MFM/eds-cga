function PosteriorOddsTable = model_comparison(ModelNames,ModelPriors,oo,options_,fname)
% Bayesian model comparison. This function computes Odds ratios and
% estimate a posterior density over a colletion of models.
%
% INPUTS
%    ModelNames       [string]     m*1 cell array of string.
%    ModelPriors      [double]     m*1 vector of prior probabilities 
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2007-2011 Dynare Team
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

NumberOfModels = size(ModelNames,2);

disp(' ')
disp(' ')
if isempty(ModelPriors)
    prior_flag = 0;% empty_prior=0
    ModelPriors = ones(NumberOfModels,1)/NumberOfModels;
else % The prior density has to sum up to one.
    prior_flag = 1;
    improper = abs(sum(ModelPriors)-1)>1e-6;
    if improper
        disp('model_comparison:: The user supplied prior distribution over models is improper...')
        disp('model_comparison:: The distribution is automatically rescaled!')
        ModelPriors=ModelPriors/sum(ModelPriors);
    end
end

% The marginal densities are based on Laplace approxiations (default) or
% modified harmonic mean estimators.
if isfield(options_,'model_comparison_approximation')
    type = options_.model_comparison_approximation;
    if strcmp(type,'Laplace')
        type = 'LaplaceApproximation';
    end
else
    type = 'LaplaceApproximation';
end

% Get the estimated logged marginal densities.
MarginalLogDensity = zeros(NumberOfModels,1);
ShortModelNames = get_short_names(ModelNames);
iname = strmatch(fname,ShortModelNames,'exact');

for i=1:NumberOfModels
    if i==iname
        mstruct.oo_ = oo;
    else
        if strcmpi(ModelNames{i}(end-3:end),'.mod') | strcmpi(ModelNames{i}(end-3:end),'.dyn')
            mstruct = load([ModelNames{i}(1:end-4) '_results.mat' ],'oo_');
        else
            mstruct = load([ModelNames{i} '_results.mat' ],'oo_');
        end
    end
    try
        eval(['MarginalLogDensity(i) = mstruct.oo_.MarginalDensity.' type ';']) 
    catch
        if strcmpi(type,'LaplaceApproximation')
            disp(['MODEL_COMPARISON: I cant''t find the Laplace approximation associated to model ' ModelNames{i}])
            return
        elseif strcmpi(type,'ModifiedHarmonicMean')
            disp(['MODEL_COMPARISON: I cant''t find the modified harmonic mean estimate associated to model ' ModelNames{i}])
            return
        end
    end
end

% In order to avoid overflow, we divide the numerator and the denominator
% of the Posterior Odds Ratio by the largest Marginal Posterior Density
lmpd = log(ModelPriors)+MarginalLogDensity;
[maxval,k] = max(lmpd);
elmpd = exp(lmpd-maxval);


% Now I display the posterior probabilities.
title = 'Model Comparison'; 
headers = char('Model',ShortModelNames{:});
if prior_flag
    labels = char('Priors','Log Marginal Density','Bayes Ratio', ...
                  'Posterior Model Probability');
    values = [ModelPriors';MarginalLogDensity';exp(lmpd-lmpd(1))'; ...
              elmpd'/sum(elmpd)];
else
    labels = char('Priors','Log Marginal Density','Bayes Ratio','Posterior Odds Ratio', ...
                  'Posterior Model Probability');
    values = [ModelPriors';MarginalLogDensity'; exp(MarginalLogDensity-MarginalLogDensity(1))'; ...
              exp(lmpd-lmpd(1))'; elmpd'/sum(elmpd)];
end

dyntable(title,headers,labels,values, 0, 15, 6);


function name = get_model_name_without_path(modelname)
idx = strfind(modelname,'\');
if isempty(idx)
    idx = strfind(modelname,'/');
end
if isempty(idx)
    name = modelname;
    return
end
name = modelname(idx(end)+1:end);


function name = get_model_name_without_extension(modelname)
idx = strfind(modelname,'.mod');
if isempty(idx)
    idx = strfind(modelname,'.dyn')
end
if isempty(idx)
    name = modelname;
    return
end
name = modelname(1:end-4);


function modellist = get_short_names(modelnames)
n = length(modelnames);
modellist = {};
for i=1:n
    name = get_model_name_without_extension(modelnames{i});
    name = get_model_name_without_path(name);
    modellist = {modellist{:} name};
end
