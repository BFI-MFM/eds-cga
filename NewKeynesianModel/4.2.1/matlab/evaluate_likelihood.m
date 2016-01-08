function [llik,parameters] = evaluate_likelihood(parameters)
% Evaluate the logged likelihood at parameters.
%
% INPUTS
%    o parameters  a string ('posterior mode','posterior mean','posterior median','prior mode','prior mean') or a vector of values for 
%                  the (estimated) parameters of the model.
%    
%    
% OUTPUTS
%    o ldens       [double]  value of the sample logged density at parameters.
%    o parameters  [double]  vector of values for the estimated parameters.
%    
% SPECIAL REQUIREMENTS
%    None
%
% REMARKS
% [1] This function cannot evaluate the likelihood of a dsge-var model...
% [2] This function use persistent variables for the dataset and the description of the missing observations. Consequently, if this function 
%     is called more than once (by changing the value of parameters) the sample *must not* change.

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

global options_ M_ bayestopt_ oo_

persistent load_data
persistent gend  data  data_index  number_of_observations  no_more_missing_observations

if nargin==0
    parameters = 'posterior mode';
end

if ischar(parameters)
    switch parameters
      case 'posterior mode'
        parameters = get_posterior_parameters('mode');
      case 'posterior mean'
        parameters = get_posterior_parameters('mean');
      case 'posterior median'
        parameters = get_posterior_parameters('median');
      case 'prior mode'
        parameters = bayestopt_.p5(:);
      case 'prior mean'
        parameters = bayestopt_.p1;
      otherwise
        disp('eval_likelihood:: If the input argument is a string, then it has to be equal to:')
        disp('                   ''posterior mode'', ')
        disp('                   ''posterior mean'', ')
        disp('                   ''posterior median'', ')
        disp('                   ''prior mode'' or')
        disp('                   ''prior mean''.')
        error
    end
end

if isempty(load_data)
    % Get the data.
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
    % Determine if a constant is needed.
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
    load_data = 1;
end

pshape_original   = bayestopt_.pshape;
bayestopt_.pshape = Inf(size(bayestopt_.pshape));
clear('priordens')%

llik = -DsgeLikelihood(parameters,gend,data,data_index,number_of_observations,no_more_missing_observations);

bayestopt_.pshape = pshape_original;