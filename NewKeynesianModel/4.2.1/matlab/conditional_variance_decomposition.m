function ConditionalVarianceDecomposition = conditional_variance_decomposition(StateSpaceModel, Steps, SubsetOfVariables,sigma_e_is_diagonal)
% This function computes the conditional variance decomposition of a given state space model
% for a subset of endogenous variables.
% 
% INPUTS 
%   StateSpaceModel     [structure]   Specification of the state space model.
%   Steps               [integer]     1*h vector of dates.
%   SubsetOfVariables   [integer]     1*q vector of indices.
%    
% OUTPUTS 
%   ConditionalVarianceDecomposition  [double] [n h p] array, where 
%                                                    n is equal to length(SubsetOfVariables)
%                                                    h is the number of Steps
%                                                    p is the number of state innovations and
% SPECIAL REQUIREMENTS
%
% [1] In this version, absence of measurement errors is assumed...

% Copyright (C) 2010 Dynare Team
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

number_of_state_innovations = ...
    StateSpaceModel.number_of_state_innovations;
transition_matrix = StateSpaceModel.transition_matrix;
number_of_state_equations = ...
    StateSpaceModel.number_of_state_equations;
order_var = StateSpaceModel.order_var;
nSteps = length(Steps);

ConditionalVariance = zeros(number_of_state_equations,nSteps,number_of_state_innovations);

if StateSpaceModel.sigma_e_is_diagonal
    B = StateSpaceModel.impulse_matrix.* ...
        repmat(sqrt(diag(StateSpaceModel.state_innovations_covariance_matrix)'),...
               number_of_state_equations,1);
else
    B = StateSpaceModel.impulse_matrix*chol(StateSpaceModel.state_innovations_covariance_matrix)';
end

for i=1:number_of_state_innovations
    BB = B(:,i)*B(:,i)';
    V = zeros(number_of_state_equations,number_of_state_equations);
    m = 1;
    for h = 1:max(Steps)
        V = transition_matrix*V*transition_matrix'+BB;
        if h == Steps(m)
            ConditionalVariance(order_var,m,i) = diag(V);
            m = m+1;
        end
    end
end

ConditionalVariance = ConditionalVariance(SubsetOfVariables,:,:);

NumberOfVariables = length(SubsetOfVariables);
SumOfVariances = zeros(NumberOfVariables,nSteps);
for h = 1:length(Steps)
    SumOfVariances(:,h) = sum(ConditionalVariance(:,h,:),3);
end

ConditionalVarianceDecomposition = zeros(NumberOfVariables,length(Steps),number_of_state_innovations); 
for i=1:number_of_state_innovations
    for h = 1:length(Steps)
        ConditionalVarianceDecomposition(:,h,i) = squeeze(ConditionalVariance(:,h,i))./SumOfVariances(:,h);
    end
end