function [xparam1, estim_params_, bayestopt_, lb, ub, M_]=set_prior(estim_params_, M_, options_)
% function [xparam1,estim_params_,bayestopt_,lb,ub]=set_prior(estim_params_)
% sets prior distributions
%
% INPUTS
%    o estim_params_    [structure] characterizing parameters to be estimated.
%    o M_               [structure] characterizing the model. 
%    o options_         [structure] 
%    
% OUTPUTS
%    o xparam1          [double]    vector of parameters to be estimated (initial values)
%    o estim_params_    [structure] characterizing parameters to be estimated
%    o bayestopt_       [structure] characterizing priors
%    o lb               [double]    vector of lower bounds for the estimated parameters. 
%    o ub               [double]    vector of upper bounds for the estimated parameters.
%    o M_               [structure] characterizing the model.
%    
% SPECIAL REQUIREMENTS
%    None

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

nvx = size(estim_params_.var_exo,1);
nvn = size(estim_params_.var_endo,1);
ncx = size(estim_params_.corrx,1);
ncn = size(estim_params_.corrn,1);
np = size(estim_params_.param_vals,1);

estim_params_.nvx = nvx;
estim_params_.nvn = nvn;
estim_params_.ncx = ncx;
estim_params_.ncn = ncn;
estim_params_.np = np;

xparam1 = [];
ub = [];
lb = [];
bayestopt_.pshape = [];
bayestopt_.p1 = []; % prior mean
bayestopt_.p2 = []; % prior standard deviation
bayestopt_.p3 = []; % lower bound
bayestopt_.p4 = []; % upper bound
bayestopt_.p5 = []; % prior mode
bayestopt_.p6 = []; % first hyper-parameter (\alpha for the BETA and GAMMA distributions, s for the INVERSE GAMMAs, expectation for the GAUSSIAN distribution, lower bound for the UNIFORM distribution).
bayestopt_.p7 = []; % second hyper-parameter (\beta for the BETA and GAMMA distributions, \nu for the INVERSE GAMMAs, standard deviation for the GAUSSIAN distribution, upper bound for the UNIFORM distribution).

bayestopt_.jscale = [];
bayestopt_.name = [];
if nvx
    xparam1 = estim_params_.var_exo(:,2);
    ub = estim_params_.var_exo(:,4); 
    lb = estim_params_.var_exo(:,3); 
    bayestopt_.pshape =  estim_params_.var_exo(:,5);
    bayestopt_.p1 =  estim_params_.var_exo(:,6);
    bayestopt_.p2 =  estim_params_.var_exo(:,7);
    bayestopt_.p3 =  estim_params_.var_exo(:,8);
    bayestopt_.p4 =  estim_params_.var_exo(:,9);
    bayestopt_.jscale =  estim_params_.var_exo(:,10);
    bayestopt_.name = cellstr(M_.exo_names(estim_params_.var_exo(:,1),:));
end
if nvn
    if isequal(M_.H,0)
        nvarobs = size(options_.varobs,1);
        M_.H = zeros(nvarobs,nvarobs);
    end
    for i=1:nvn
        obsi_ = strmatch(deblank(M_.endo_names(estim_params_.var_endo(i,1),:)),deblank(options_.varobs),'exact');
        if isempty(obsi_)
            error(['The variable ' deblank(M_.endo_names(estim_params_.var_endo(i,1),:)) ' has to be declared as observable since you assume a measurement error on it.'])
        end
        estim_params_.var_endo(i,1) = obsi_;
    end
    xparam1 = [xparam1; estim_params_.var_endo(:,2)];
    ub = [ub; estim_params_.var_endo(:,4)]; 
    lb = [lb; estim_params_.var_endo(:,3)]; 
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.var_endo(:,5)];
    bayestopt_.p1 = [ bayestopt_.p1; estim_params_.var_endo(:,6)];
    bayestopt_.p2 = [ bayestopt_.p2; estim_params_.var_endo(:,7)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.var_endo(:,8)];
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.var_endo(:,9)];
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.var_endo(:,10)];
    if isempty(bayestopt_.name)
        bayestopt_.name = cellstr(char(options_.varobs(estim_params_.var_endo(:,1),:)));
    else
        bayestopt_.name = cellstr(char(char(bayestopt_.name), options_.varobs(estim_params_.var_endo(:,1),:)));
    end
end
if ncx
    xparam1 = [xparam1; estim_params_.corrx(:,3)];
    ub = [ub; max(min(estim_params_.corrx(:,5),1),-1)];
    lb = [lb; max(min(estim_params_.corrx(:,4),1),-1)];
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.corrx(:,6)];
    bayestopt_.p1 = [ bayestopt_.p1; estim_params_.corrx(:,7)];
    bayestopt_.p2 = [ bayestopt_.p2; estim_params_.corrx(:,8)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.corrx(:,9)];
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.corrx(:,10)];
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.corrx(:,11)];
    if isempty(bayestopt_.name)
        bayestopt_.name = cellstr(char(char(strcat(cellstr(M_.exo_names(estim_params_.corrx(:,1),:)), ...
                                                   ',' , cellstr(M_.exo_names(estim_params_.corrx(:,2),:))))));
    else
        bayestopt_.name = cellstr(char(char(bayestopt_.name), char(strcat(cellstr(M_.exo_names(estim_params_.corrx(:,1),:)), ...
                                                          ',' , cellstr(M_.exo_names(estim_params_.corrx(:,2),:))))));
    end
end
if ncn
    if isequal(M_.H,0)
        nvarobs = size(options_.varobs,1);
        M_.H = zeros(nvarobs,nvarobs);
    end
    xparam1 = [xparam1; estim_params_.corrn(:,3)];
    ub = [ub; max(min(estim_params_.corrn(:,5),1),-1)];
    lb = [lb; max(min(estim_params_.corrn(:,4),1),-1)];
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.corrn(:,6)];
    bayestopt_.p1 = [ bayestopt_.p1; estim_params_.corrn(:,7)];
    bayestopt_.p2 = [ bayestopt_.p2; estim_params_.corrn(:,8)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.corrn(:,9)];
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.corrn(:,10)];
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.corrn(:,11)];
    if isempty(bayestopt_.name)
        bayestopt_.name = cellstr(char(char(strcat(cellstr(M_.endo_names(estim_params_.corrn(:,1),:)),...
                                                   ',' ,  cellstr(M_.endo_names(estim_params_.corrn(:,2),:))))));
    else
        bayestopt_.name = cellstr(char(char(bayestopt_.name), char(strcat(cellstr(M_.endo_names(estim_params_.corrn(:,1),:)),...
                                                          ',' ,  cellstr(M_.endo_names(estim_params_.corrn(:,2),:))))));
    end
end
if np
    xparam1 = [xparam1; estim_params_.param_vals(:,2)];
    ub = [ub; estim_params_.param_vals(:,4)];
    lb = [lb; estim_params_.param_vals(:,3)];
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.param_vals(:,5)];
    bayestopt_.p1 = [ bayestopt_.p1; estim_params_.param_vals(:,6)];
    bayestopt_.p2 = [ bayestopt_.p2; estim_params_.param_vals(:,7)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.param_vals(:,8)];
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.param_vals(:,9)];
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.param_vals(:,10)];
    if isempty(bayestopt_.name)
        bayestopt_.name = cellstr(char(M_.param_names(estim_params_.param_vals(:,1),:)));
    else
        bayestopt_.name = cellstr(char(char(bayestopt_.name),M_.param_names(estim_params_.param_vals(:,1),:)));
    end
end

bayestopt_.ub = ub;
bayestopt_.lb = lb;

bayestopt_.p6 = NaN(size(bayestopt_.p1)) ;
bayestopt_.p7 = bayestopt_.p6 ;

% generalized location parameters by default for beta distribution
k = find(bayestopt_.pshape == 1);
k1 = find(isnan(bayestopt_.p3(k)));
bayestopt_.p3(k(k1)) = zeros(length(k1),1);
k1 = find(isnan(bayestopt_.p4(k)));
bayestopt_.p4(k(k1)) = ones(length(k1),1);
for i=1:length(k)
    if (bayestopt_.p1(k(i))<bayestopt_.p3(k(i))) || (bayestopt_.p1(k(i))>bayestopt_.p4(k(i)))
        error(['The prior mean of ' bayestopt_.name{k(i)} ' has to be between the lower (' num2str(bayestopt_.p3(k(i))) ') and upper (' num2str(bayestopt_.p4(k(i))) ') bounds of the beta prior density!']);
    end
    mu = (bayestopt_.p1(k(i))-bayestopt_.p3(k(i)))/(bayestopt_.p4(k(i))-bayestopt_.p3(k(i)));
    stdd = bayestopt_.p2(k(i))/(bayestopt_.p4(k(i))-bayestopt_.p3(k(i)));
    bayestopt_.p6(k(i)) = (1-mu)*mu^2/stdd^2 - mu ;
    bayestopt_.p7(k(i)) = bayestopt_.p6(k(i))*(1/mu-1) ;
    m = compute_prior_mode([ bayestopt_.p6(k(i)) , bayestopt_.p7(k(i)) , bayestopt_.p3(k(i)) , bayestopt_.p4(k(i)) ],1);
    if length(m)==1
        bayestopt_.p5(k(i)) = m;
    else
        disp(['Prior distribution for parameter ' bayestopt_.name{k(i)}  ' has two modes!'])
        bayestopt_.p5(k(i)) = bayestopt_.p1(k(i)) ; 
    end
end

% generalized location parameter by default for gamma distribution
k =  find(bayestopt_.pshape == 2);
k1 = find(isnan(bayestopt_.p3(k)));
k2 = find(isnan(bayestopt_.p4(k)));
bayestopt_.p3(k(k1)) = zeros(length(k1),1);
bayestopt_.p4(k(k2)) = Inf(length(k2),1);
for i=1:length(k)
    if isinf(bayestopt_.p2(k(i)))
        error(['Infinite prior standard deviation for parameter ' bayestopt_.name{k(i)}  ' is not allowed (Gamma prior)!'])
    end
    mu = bayestopt_.p1(k(i))-bayestopt_.p3(k(i));
    bayestopt_.p7(k(i)) = bayestopt_.p2(k(i))^2/mu ;
    bayestopt_.p6(k(i)) = mu/bayestopt_.p7(k(i)) ;  
    bayestopt_.p5(k(i)) = compute_prior_mode([ bayestopt_.p6(k(i)) , bayestopt_.p7(k(i)) , bayestopt_.p3(k(i)) ], 2) ;
end

% truncation parameters by default for normal distribution
k  = find(bayestopt_.pshape == 3);
k1 = find(isnan(bayestopt_.p3(k)));
bayestopt_.p3(k(k1)) = -Inf*ones(length(k1),1);
k1 = find(isnan(bayestopt_.p4(k)));
bayestopt_.p4(k(k1)) = Inf*ones(length(k1),1);
for i=1:length(k)
    bayestopt_.p6(k(i)) = bayestopt_.p1(k(i)) ; 
    bayestopt_.p7(k(i)) = bayestopt_.p2(k(i)) ;
    bayestopt_.p5(k(i)) = bayestopt_.p1(k(i)) ;
end

% inverse gamma distribution (type 1)
k = find(bayestopt_.pshape == 4);
k1 = find(isnan(bayestopt_.p3(k)));
k2 = find(isnan(bayestopt_.p4(k)));
bayestopt_.p3(k(k1)) = zeros(length(k1),1);
bayestopt_.p4(k(k2)) = Inf(length(k2),1);
for i=1:length(k)
    [bayestopt_.p6(k(i)),bayestopt_.p7(k(i))] = ...
        inverse_gamma_specification(bayestopt_.p1(k(i))-bayestopt_.p3(k(i)),bayestopt_.p2(k(i)),1) ;
    bayestopt_.p5(k(i)) = compute_prior_mode([ bayestopt_.p6(k(i)) , bayestopt_.p7(k(i)) , bayestopt_.p3(k(i)) ], 4) ;
end  

% uniform distribution
k = find(bayestopt_.pshape == 5);
for i=1:length(k)
    [bayestopt_.p1(k(i)),bayestopt_.p2(k(i)),bayestopt_.p6(k(i)),bayestopt_.p7(k(i))] = ...
        uniform_specification(bayestopt_.p1(k(i)),bayestopt_.p2(k(i)),bayestopt_.p3(k(i)),bayestopt_.p4(k(i)));
    bayestopt_.p3(k(i)) = bayestopt_.p6(k(i)) ;
    bayestopt_.p4(k(i)) = bayestopt_.p7(k(i)) ;
    bayestopt_.p5(k(i)) = NaN ;
end

% inverse gamma distribution (type 2)
k = find(bayestopt_.pshape == 6);
k1 = find(isnan(bayestopt_.p3(k)));
k2 = find(isnan(bayestopt_.p4(k)));
bayestopt_.p3(k(k1)) = zeros(length(k1),1);
bayestopt_.p4(k(k2)) = Inf(length(k2),1);
for i=1:length(k)
    [bayestopt_.p6(k(i)),bayestopt_.p7(k(i))] = ...
        inverse_gamma_specification(bayestopt_.p1(k(i))-bayestopt_.p3(k(i)),bayestopt_.p2(k(i)),2);
    bayestopt_.p5(k(i)) = compute_prior_mode([ bayestopt_.p6(k(i)) , bayestopt_.p7(k(i)) , bayestopt_.p3(k(i)) ], 6) ;
end

k = find(isnan(xparam1));
xparam1(k) = bayestopt_.p1(k);

% I create subfolder M_.dname/prior if needed.
CheckPath('prior');

% I save the prior definition if the prior has changed.
if exist([ M_.dname '/prior/definition.mat'])
    old = load([M_.dname '/prior/definition.mat'],'bayestopt_');
    prior_has_changed = 0;
    if length(bayestopt_.p1)==length(old.bayestopt_.p1)
        if any(bayestopt_.p1-old.bayestopt_.p1)
            prior_has_changed = 1;
        end
        if any(bayestopt_.p2-old.bayestopt_.p2)
            prior_has_changed = 1;
        end
        if any(bayestopt_.p3-old.bayestopt_.p3)
            prior_has_changed = 1;
        end
        if any(bayestopt_.p4-old.bayestopt_.p4)
            prior_has_changed = 1;
        end
        if any(bayestopt_.p5-old.bayestopt_.p5)
            prior_has_changed = 1;
        end
        if any(bayestopt_.p6-old.bayestopt_.p6)
            prior_has_changed = 1;
        end
        if any(bayestopt_.p7-old.bayestopt_.p7)
            prior_has_changed = 1;
        end
    else
        prior_has_changed = 1;
    end
    if prior_has_changed
        delete([M_.dname '/prior/definition.mat']);
        save([M_.dname '/prior/definition.mat'],'bayestopt_');
    end
else
    save([M_.dname '/prior/definition.mat'],'bayestopt_');
end

% initialize persistent variables in priordens()
priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7, ...
          bayestopt_.p3,bayestopt_.p4,1);

% Put bayestopt_ in matlab's workspace
assignin('base','bayestopt_',bayestopt_);