function oo_ = compute_moments_varendo(type,options_,M_,oo_,var_list_)
% Computes the second order moments (autocorrelation function, covariance
% matrix and variance decomposition) distributions for all the endogenous variables selected in
% var_list_. The results are saved in oo_
%  
% INPUTS:
%   type            [string]       'posterior' or 'prior'
%   options_        [structure]    Dynare structure.
%   M_              [structure]    Dynare structure (related to model definition).
%   oo_             [structure]    Dynare structure (results).
%   var_list_       [string]       Array of string with endogenous variable names.
%    
% OUTPUTS
%   oo_             [structure]    Dynare structure (results).
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2008-2010 Dynare Team
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

if strcmpi(type,'posterior')
    posterior = 1;
    if nargin==4
        var_list_ = options_.varobs;
    end
elseif strcmpi(type,'prior')
    posterior = 0;
    if nargin==4
        var_list_ = options_.prior_analysis_endo_var_list;
        if isempty(var_list_)
            options_.prior_analysis_var_list = options_.varobs;
        end
    end
else
    disp('compute_moments_varendo:: Unknown type!')
    error()
end

NumberOfEndogenousVariables = rows(var_list_);
NumberOfExogenousVariables = M_.exo_nbr;
list_of_exogenous_variables = M_.exo_names;
NumberOfLags = options_.ar;
if isfield(options_,'conditional_variance_decomposition')
    Steps = options_.conditional_variance_decomposition;
else
    Steps = 0;
end

% COVARIANCE MATRIX.
if posterior
    for i=1:NumberOfEndogenousVariables
        for j=i:NumberOfEndogenousVariables
            oo_ = posterior_analysis('variance',var_list_(i,:),var_list_(j,:),[],options_,M_,oo_);
        end
    end
else
    for i=1:NumberOfEndogenousVariables
        for j=i:NumberOfEndogenousVariables
            oo_ = prior_analysis('variance',var_list_(i,:),var_list_(j,:),[],options_,M_,oo_);
        end
    end
end
% CORRELATION FUNCTION.
if posterior
    for h=NumberOfLags:-1:1
        for i=1:NumberOfEndogenousVariables
            for j=1:NumberOfEndogenousVariables
                oo_ = posterior_analysis('correlation',var_list_(i,:),var_list_(j,:),h,options_,M_,oo_);
            end
        end
    end
else
    for h=NumberOfLags:-1:1
        for i=1:NumberOfEndogenousVariables
            for j=1:NumberOfEndogenousVariables
                oo_ = prior_analysis('correlation',var_list_(i,:),var_list_(j,:),h,options_,M_,oo_);
            end
        end
    end
end
% VARIANCE DECOMPOSITION.
if M_.exo_nbr > 1
    if posterior
        for i=1:NumberOfEndogenousVariables
            for j=1:NumberOfExogenousVariables
                oo_ = posterior_analysis('decomposition',var_list_(i,:),M_.exo_names(j,:),[],options_,M_,oo_);
            end
        end
    else
        for i=1:NumberOfEndogenousVariables
            for j=1:NumberOfExogenousVariables
                oo_ = prior_analysis('decomposition',var_list_(i,:),M_.exo_names(j,:),[],options_,M_,oo_);
            end
        end        
    end
    % CONDITIONAL VARIANCE DECOMPOSITION.
    if Steps
        if posterior
            for i=1:NumberOfEndogenousVariables
                for j=1:NumberOfExogenousVariables
                    oo_ = posterior_analysis('conditional decomposition',i,M_.exo_names(j,:),Steps,options_,M_,oo_);
                end
            end
        else
            for i=1:NumberOfEndogenousVariables
                for j=1:NumberOfExogenousVariables
                    oo_ = prior_analysis('conditional decomposition',var_list_(i,:),M_.exo_names(j,:),Steps,options_,M_,oo_);
                end
            end
        end
    end
end
