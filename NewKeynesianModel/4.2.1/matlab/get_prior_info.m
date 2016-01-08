function get_prior_info(info)
% Computes various prior statistics.
%  
% INPUTS
%   info     [integer]   scalar specifying what has to be done.
%  
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2009 Dynare Team
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
global options_ M_ estim_params_ oo_

if ~nargin
    info = 0;
end

M_.dname = M_.fname;

[xparam1,estim_params_,bayestopt_,lb,ub,M_] = set_prior(estim_params_,M_,options_);
plot_priors(bayestopt_,M_,options_);

PriorNames = { 'Beta' , 'Gamma' , 'Gaussian' , 'Inverted Gamma' , 'Uniform' , 'Inverted Gamma -- 2' };

if size(M_.param_names,1)==size(M_.param_names_tex,1)% All the parameters have a TeX name.
    fidTeX = fopen('priors_data.tex','w+');
    % Column 1: a string for the name of the prior distribution. 
    % Column 2: the prior mean.
    % Column 3: the prior standard deviation.
    % Column 4: the lower bound of the prior density support.
    % Column 5: the upper bound of the prior density support.
    % Column 6: the lower bound of the interval containing 80% of the prior mass. 
    % Column 7: the upper bound of the interval containing 80% of the prior mass.
    prior_trunc_backup = options_.prior_trunc ;
    options_.prior_trunc = (1-options_.prior_interval)/2 ;
    PriorIntervals = prior_bounds(bayestopt_) ;
    options_.prior_trunc = prior_trunc_backup ;
    for i=1:size(bayestopt_.name,1)
        [tmp,TexName] = get_the_name(i,1);
        PriorShape = PriorNames{ bayestopt_.pshape(i) };
        PriorMean = bayestopt_.p1(i);
        PriorStandardDeviation = bayestopt_.p2(i);
        switch bayestopt_.pshape(i)
          case { 1 , 5 }
            LowerBound = bayestopt_.p3(i);
            UpperBound = bayestopt_.p4(i);
          case { 2 , 4 , 6 }
            LowerBound = bayestopt_.p3(i);
            UpperBound = '$\infty$';
          case 3
            if isinf(bayestopt_.p3(i))
                LowerBound = '$-\infty$';
            else
                LowerBound = bayestopt_.p3(i);
            end
            if isinf(bayestopt_.p4(i))
                UpperBound = '$\infty$';
            else
                UpperBound = bayestopt_.p4(i);
            end
          otherwise
            error('get_prior_info:: Dynare bug!')
        end
        format_string = build_format_string(bayestopt_,i);
        fprintf(fidTeX,format_string, ...
                TexName, ...
                PriorShape, ...
                PriorMean, ...
                PriorStandardDeviation, ...
                LowerBound, ...
                UpperBound, ...
                PriorIntervals(i,1), ...
                PriorIntervals(i,2) );
    end
    fclose(fidTeX);
end

M_.dname = M_.fname;

if info==1% Prior simulations (BK).
    results = prior_sampler(0,M_,bayestopt_,options_,oo_);
    % Display prior mass info
    disp(['Prior mass = ' num2str(results.prior.mass)])
    disp(['BK indeterminacy share                = ' num2str(results.bk.indeterminacy_share)])
    disp(['BK unstability share                  = ' num2str(results.bk.unstability_share)])
    disp(['BK singularity share                  = ' num2str(results.bk.singularity_share)])
    disp(['Complex jacobian share                = ' num2str(results.jacobian.problem_share)])
    disp(['mjdgges crash share                   = ' num2str(results.dll.problem_share)])
    disp(['Steady state problem share            = ' num2str(results.ss.problem_share)])
    disp(['Complex steady state  share           = ' num2str(results.ss.complex_share)])
    disp(['Analytical steady state problem share = ' num2str(results.ass.problem_share)])
end

if info==2% Prior optimization.
          % Initialize to the prior mode if possible
    k = find(~isnan(bayestopt_.p5));
    xparam1(k) = bayestopt_.p5(k);
    % Pertubation of the initial condition.
    look_for_admissible_initial_condition = 1;
    scale = 1.0;
    iter  = 0;
    while look_for_admissible_initial_condition
        xinit = xparam1+scale*randn(size(xparam1));
        if all(xinit>bayestopt_.p3) && all(xinit<bayestopt_.p4)
            look_for_admissible_initial_condition = 0;
        else
            if iter == 2000;
                scale = scale/1.1;
                iter  = 0;
            else
                iter = iter+1;
            end
        end
    end
    % Maximization
    [xparams,lpd,hessian] = ...
        maximize_prior_density(xinit, bayestopt_.pshape, ...
                               bayestopt_.p6, ...
                               bayestopt_.p7, ...
                               bayestopt_.p3, ...
                               bayestopt_.p4);
    % Display the results.
    disp(' ')
    disp(' ')
    disp('------------------')
    disp('PRIOR OPTIMIZATION')
    disp('------------------')
    disp(' ')
    for i = 1:length(xparams)
        disp(['deep parameter ' int2str(i) ': ' get_the_name(i,0) '.'])
        disp(['  Initial condition ....... ' num2str(xinit(i)) '.'])
        disp(['  Prior mode .............. ' num2str(bayestopt_.p5(i)) '.'])
        disp(['  Optimized prior mode .... ' num2str(xparams(i)) '.'])
        disp(' ')
    end
end

if info==3% Prior simulations (2nd order moments).
    oo_ = compute_moments_varendo('prior',options_,M_,oo_);
end 


function format_string = build_format_string(bayestopt,i)
format_string = ['%s & %s & %6.4f &'];
if isinf(bayestopt.p2(i))
    format_string = [ format_string , ' %s &'];
else
    format_string = [ format_string , ' %6.4f &'];
end
if isinf(bayestopt.p3(i))
    format_string = [ format_string , ' %s &'];
else
    format_string = [ format_string , ' %6.4f &'];
end
if isinf(bayestopt.p4(i))
    format_string = [ format_string , ' %s &'];
else
    format_string = [ format_string , ' %6.4f &'];
end
format_string = [ format_string , ' %6.4f & %6.4f \\\\ \n'];