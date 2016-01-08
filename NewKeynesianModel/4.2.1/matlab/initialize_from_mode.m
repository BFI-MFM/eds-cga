function estim_params_ = initialize_from_mode(fname,M_,estim_params_)
% function estim_params_ = initialize_from_mode(fname,M_,estim_params_)
% initialize parameters and initial value of estimated parameters
% from a *_mode.mat file    
%  
% INPUTS
%   fname:  mode file name (*.mat file)
%   M_:     sructure of model characteristics
%   estim_params_: structure of estimated parameters
%  
% OUTPUTS
%   estim_params:  modified structure of estimated parameters
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2003-2009 Dynare Team
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


load(fname,'xparam1','parameter_names');

endo_names = M_.endo_names;
exo_names = M_.exo_names;
param_names = M_.param_names;
param_vals = estim_params_.param_vals;
var_exo = estim_params_.var_exo;
var_endo = estim_params_.var_endo;
corrx = estim_params_.corrx;
corrn = estim_params_.corrn;
for i=1:length(parameter_names)
    name = parameter_names{i};
    k1 = strmatch(name,param_names,'exact');
    if ~isempty(k1)
        k2 = find(param_vals(:,1) == k1);
        if ~isempty(k2)
            estim_params_.param_vals(k2,2) = xparam1(i);
        end
        M_.params(i) = xparam1(i);
        continue
    end
    k3 = strfind(name,',');
    if isempty(k3)
        k1 = strmatch(name,exo_names,'exact');
        if ~isempty(k1)
            k2 = find(var_exo(:,1) == k1);
            if ~isempty(k2)
                estim_params_.var_exo(k2,2) = xparam1(i);
            end
            M_.Sigma_e(k1,k1) = xparam1(i)^2;
            continue
        end
        k1 = strmatch(name,endo_names,'exact');
        if ~isempty(k1)
            k2 = find(var_endo(:,1) == k1);
            if ~isempty(k2)
                estim_params_.var_endo(k2,2) = xparam1(i);
            end
            M_.H(k1,k1) = xparam1(i)^2;
            continue
        end
    else
        k1 = strmatch(name(1:k3-1),exo_names,'exact');
        k1a = strmatch(name(k3+1:end),exo_names,'exact');
        if ~isempty(k1) & ~isempty(k1a)
            k2 = find(corrx(:,1) == k1 & corrx(:,2) == k1a);
            if ~isempty(k2)
                estim_params_.corrx(k2,3) = xparam1(i);
            end
            M_.Sigma_e(k1,k1a) = xparam1(i)*sqrt(M_.Sigma_e(k1,k1)+M_.Sigma_e(k1a,k1a));
            M_.Sigma_e(k1a,k1) = M_.Sigma_e(k1,k1a);
            continue
        end
        k1 = strmatch(name(1:k3-1),endo_names,'exact');
        k1a = strmatch(name(k3+1:end),endo_names,'exact');
        if ~isempty(k1) & ~isempty(k1a)
            k2 = find(corrn(:,1) == k1 & corrn(:,2) == k1a);
            if ~isempty(k2)
                estim_params_.corrn(k2,3) = xparam1(i);
            end
            M_.H(k1,k1a) = xparam1(i)*sqrt(M_.H(k1,k1)+M_.H(k1a,k1a));
            M_.H(k1a,k1) = M_.H(k1,k1a);
            continue
        end
    end
    error([name 'doesn''t exist in this model'])
end

