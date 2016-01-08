function initvalf(fname_)
% function initvalf(fname_)
%
% Reads an initial path from the 'fname_' file for exogenous and endogenous variables   
%
% INPUTS
%    fname_:         name of the function or file containing the data
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    All variables local to this function have an underscore appended to
%    their name, to minimize clashes with model variables loaded by this function.

% Copyright (C) 2003-2010 Dynare Team
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

global M_ oo_ options_

series_ = 1;
if exist(fname_) == 2
    eval(fname_);
elseif exist([fname_ '.xls']) == 2
    [data_,names_v_]=xlsread([fname_ '.xls']);
    series_ = 0;
elseif exist([fname_ '.mat']) == 2
    load(fname_);
end

options_.initval_file = 1;
oo_.endo_simul = [];
oo_.exo_simul = [];

for i_=1:size(M_.endo_names,1)
    if series_ == 1
        x_ = eval(M_.endo_names(i_,:));
        oo_.endo_simul = [oo_.endo_simul; x_'];
    else
        k_ = strmatch(upper(M_.endo_names(i_,:)),names_v_,'exact');
        if isempty(k_)
            error(['INITVAL_FILE: ' M_.endo_names(i_,:) ' not found'])
        end
        x_ = data_(:,k_);
        oo_.endo_simul = [oo_.endo_simul; x_']; 
    end
end

for i_=1:size(M_.exo_names,1)
    if series_ == 1
        x_ = eval(M_.exo_names(i_,:) );
        oo_.exo_simul = [oo_.exo_simul x_];
    else
        k_ = strmatch(upper(M_.exo_names(i_,:)),names_v_,'exact');
        if isempty(k_)
            error(['INITVAL_FILE: ' M_.exo_names(i_,:) ' not found'])
        end
        x_ = data_(:,k_);
        oo_.exo_simul = [oo_.exo_simul x_]; 
    end
end
