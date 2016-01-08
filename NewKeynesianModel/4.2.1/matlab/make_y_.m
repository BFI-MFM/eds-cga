function make_y_
% function make_y_
% forms oo_.endo_simul as guess values for deterministic simulations
%  
% INPUTS
%   ...
% OUTPUTS
%   ...
% ALGORITHM
%   ...
% SPECIAL REQUIREMENTS
%   none
%  

% Copyright (C) 1996-2009 Dynare Team
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

global M_ options_ oo_ ys0_ 

options_ = set_default_option(options_,'periods',0);

if isempty(oo_.steady_state)
    oo_.steady_state = zeros(M_.endo_nbr,1);
end

if isempty(oo_.endo_simul)
    if isempty(ys0_)
        oo_.endo_simul = [oo_.steady_state*ones(1,M_.maximum_lag+options_.periods+M_.maximum_lead)];
    else
        oo_.endo_simul = [ys0_*ones(1,M_.maximum_lag) oo_.steady_state*ones(1,options_.periods+M_.maximum_lead)];
    end
elseif size(oo_.endo_simul,2) < M_.maximum_lag+M_.maximum_lead+options_.periods
    switch options_.deterministic_simulation_initialization 
      case 0
        oo_.endo_simul = [oo_.endo_simul ...
                          oo_.steady_state*ones(1,M_.maximum_lag+options_.periods+M_.maximum_lead-size(oo_.endo_simul,2),1)];
      case 1% A linear approximation is used to initialize the solution.
        oldopt = options_;
        options_.order = 1;
        dr = oo_.dr;
        dr.ys = oo_.steady_state;
        [dr,info,M_,options_,oo_]=dr1(dr,0,M_,options_,oo_);
        exogenous_variables = zeros(M_.maximum_lag+options_.periods+M_.maximum_lead-size(oo_.endo_simul,2)+1,0);
        y0 = oo_.endo_simul(:,1:M_.maximum_lag);
        oo_.endo_simul=simult_(y0,dr,exogenous_variables,1);
        options_ = oldopt;
      case 2% Homotopic mod: Leave endo_simul as it is.
      otherwise
        error('Unknown method.')
    end
end