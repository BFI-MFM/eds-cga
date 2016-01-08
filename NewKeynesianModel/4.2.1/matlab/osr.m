function osr(var_list,params,i_var,W)

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

global M_ options_ oo_  

options_.order = 1;
options_ = set_default_option(options_,'linear',0);
options_ = set_default_option(options_,'ar',5);
options_ = set_default_option(options_,'irf',40);
options_ = set_default_option(options_,'drop',100);
options_ = set_default_option(options_,'replic',1);
options_ = set_default_option(options_,'nomoments',0);
options_ = set_default_option(options_,'nocorr',0);
options_ = set_default_option(options_,'hp_filter',0);
options_ = set_default_option(options_,'hp_ngrid',512);
options_ = set_default_option(options_,'simul',0);
options_ = set_default_option(options_,'periods',1);

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

make_ex_;

np = size(params,1);
i_params = zeros(np,1);
for i=1:np
    i_params(i) = strmatch(deblank(params(i,:)),M_.param_names,'exact');
end

disp(' ')
disp('OPTIMAL SIMPLE RULE')
disp(' ')
osr1(i_params,i_var,W);

stoch_simul(var_list);