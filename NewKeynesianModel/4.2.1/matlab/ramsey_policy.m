function info = ramsey_policy(var_list)

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

global options_ oo_ M_

options_.ramsey_policy = 1;
oldoptions = options_;
options_.order = 1;
info = stoch_simul(var_list);

if options_.noprint == 0
    disp_steady_state(M_,oo_)
    for i=M_.orig_endo_nbr:M_.endo_nbr
        if strmatch('mult_',M_.endo_names(i,:))
            disp(sprintf('%s \t\t %g',M_.endo_names(i,:), ...
                         oo_.dr.ys(i)));
        end
    end
end


%oo_ = evaluate_planner_objective(oo_.dr,M_,oo_,options_);

options_ = oldoptions;