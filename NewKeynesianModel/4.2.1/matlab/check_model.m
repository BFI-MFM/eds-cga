function check_model()

% Copyright (C) 2005-2011 Dynare Team
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

global M_

xlen = M_.maximum_exo_lag+M_.maximum_exo_lead + 1;
if ~ M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0
    error ('RESOL: Error in model specification: some variables don"t appear as current') ;
end

if xlen > 1
    error (['RESOL: stochastic exogenous variables must appear only at the' ...
            ' current period. Use additional endogenous variables']) ;
end

if (M_.exo_det_nbr > 0) && (M_.maximum_lag > 1 || M_.maximum_lead > 1)
    error(['Exogenous deterministic variables are currently only allowed in' ...
           ' models with leads and lags on only one period'])
end

