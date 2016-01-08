function make_ex_
% forms oo_.exo_simul and oo_.exo_det_simul
%
% INPUTS
%   none
%    
% OUTPUTS
%   none
%
% ALGORITHM
%   
% SPECIAL REQUIREMENTS
%  

% Copyright (C) 1996-2011 Dynare Team
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

global M_ options_ oo_ ex0_ ex_det0_

options_ = set_default_option(options_,'periods',0);

if isempty(oo_.exo_steady_state)
    oo_.exo_steady_state = zeros(M_.exo_nbr,1);
end
if M_.exo_det_nbr > 1 && isempty(oo_.exo_det_steady_state)
    oo_.exo_det_steady_state = zeros(M_.exo_det_nbr,1);
end
if isempty(oo_.exo_simul)
    if isempty(ex0_)
        oo_.exo_simul = repmat(oo_.exo_steady_state',M_.maximum_lag+options_.periods+M_.maximum_lead,1);
    else
        oo_.exo_simul = [ repmat(ex0_',M_.maximum_lag,1) ; repmat(oo_.exo_steady_state',options_.periods+M_.maximum_lead,1) ];
    end
elseif size(oo_.exo_simul,1) < M_.maximum_lag+M_.maximum_lead+options_.periods
    oo_.exo_simul = [ oo_.exo_simul ; repmat(oo_.exo_steady_state',M_.maximum_lag+options_.periods+M_.maximum_lead-size(oo_.exo_simul,1),1) ];
end

if M_.exo_det_nbr > 0
    if isempty(oo_.exo_det_simul)
        if isempty(ex_det0_)
            oo_.exo_det_simul = [ones(M_.maximum_lag+options_.periods+M_.maximum_lead,1)*oo_.exo_det_steady_state'];
        else
            oo_.exo_det_simul = [ones(M_.maximum_lag,1)*ex_det0_';ones(options_.periods+M_.maximum_lead,1)*oo_.exo_det_steady_state'];
        end
    elseif size(oo_.exo_det_simul,1) < M_.maximum_lag+M_.maximum_lead+options_.periods
        oo_.exo_det_simul = [oo_.exo_det_simul; ones(M_.maximum_lag+options_.periods+M_.maximum_lead-size(oo_.exo_det_simul,1),1)*oo_.exo_det_steady_state'];
    end
end