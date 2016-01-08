function [cmm, mm] = simulated_moment_uncertainty(indx, periods, replic)

%
% Copyright (C) 2009-2010 Dynare Team
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

global options_ oo_

mm=zeros(length(indx),replic);

disp('Evaluting simulated moment uncertainty ... please wait')
disp(['Doing ',int2str(replic),' replicas of length ',int2str(periods),' periods.'])
noprint0 = options_.noprint;
for j=1:replic;
    options_.irf = 0;
    options_.noprint = 1;
    options_.order = 1;
    options_.periods = periods;
    info = stoch_simul(options_.varobs);
    dum=[oo_.mean; dyn_vech(oo_.var)];
    sd = sqrt(diag(oo_.var));
    for i=1:options_.ar;
        dum=[dum; vec(oo_.autocorr{i}.*(sd*sd'))];
    end
    mm(:,j)=dum(indx);
end;

options_.noprint = noprint0;
cmm = cov(mm');
disp('Simulated moment uncertainty ... done!')
