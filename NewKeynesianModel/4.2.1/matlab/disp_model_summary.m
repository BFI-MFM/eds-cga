function disp_model_summary(M,dr)

% function disp_model_summary(M)
%     displays the model summary
%
% INPUTS
%   M         [matlab structure] Definition of the model.           
%   dr        [matlab structure] Decision rules
%
% Copyright (C) 2001-2010 Dynare Team
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

disp(' ')
disp('MODEL SUMMARY')
disp(' ')
disp(['  Number of variables:         ' int2str(M.endo_nbr)])
disp(['  Number of stochastic shocks: ' int2str(M.exo_nbr)])
disp(['  Number of state variables:   ' ...
      int2str(length(find(dr.kstate(:,2) <= M.maximum_lag+1)))])
disp(['  Number of jumpers:           ' ...
      int2str(length(find(dr.kstate(:,2) == M.maximum_lag+2)))])
disp(['  Number of static variables:  ' int2str(dr.nstatic)])
my_title='MATRIX OF COVARIANCE OF EXOGENOUS SHOCKS';
labels = deblank(M.exo_names);
headers = char('Variables',labels);
lh = size(labels,2)+2;
dyntable(my_title,headers,labels,M.Sigma_e,lh,10,6);

