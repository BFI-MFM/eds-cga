function steady()
% function steady()
% computes and prints the steady state calculations
%  
% INPUTS
%   none
%  
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

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

global M_ oo_ options_ ys0_ 

test_for_deep_parameters_calibration(M_);

options_ = set_default_option(options_,'jacobian_flag',1);
options_ = set_default_option(options_,'steadystate_flag',0);
if exist([M_.fname '_steadystate.m'])
    options_.steadystate_flag = 1;
end 

if options_.steadystate_flag && options_.homotopy_mode
    error('STEADY: Can''t use homotopy when providing a steady state external file');
end

switch options_.homotopy_mode
  case 1
    homotopy1(options_.homotopy_values, options_.homotopy_steps);
  case 2
    homotopy2(options_.homotopy_values, options_.homotopy_steps);
  case 3
    homotopy3(options_.homotopy_values, options_.homotopy_steps);
end

steady_;

disp_steady_state(M_,oo_);

if isempty(ys0_)
    oo_.endo_simul(:,1:M_.maximum_lag) = oo_.steady_state * ones(1, M_.maximum_lag);
    %%% Unless I'm wrong, this is (should be?) done in make_y_.m
% $$$   else
% $$$     options_ =set_default_option(options_,'periods',1);
% $$$     oo_.endo_simul(:,M_.maximum_lag+1:M_.maximum_lag+options_.periods+ ...
% $$$                    M_.maximum_lead) = oo_.steady_state * ones(1,options_.periods+M_.maximum_lead);
end