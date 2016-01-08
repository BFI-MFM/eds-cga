function [posterior_mean,posterior_covariance,posterior_mode,posterior_kernel_at_the_mode] = compute_mh_covariance_matrix()
% Estimation of the posterior covariance matrix, posterior mean, posterior mode and evaluation of the posterior kernel at the
% estimated mode, using the draws from a metropolis-hastings. The estimated posterior mode and covariance matrix are saved in
% a file <M_.fname>_mh_mode.mat. 
% 
% INPUTS 
%   None.
%  
% OUTPUTS
%   o  posterior_mean                [double]  (n*1) vector, posterior expectation of the parameters.
%   o  posterior_covariance          [double]  (n*n) matrix, posterior covariance of the parameters (computed from previous metropolis hastings).
%   o  posterior_mode                [double]  (n*1) vector, posterior mode of the parameters. 
%   o  posterior_kernel_at_the_mode  [double]  scalar.
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2006-2010 Dynare Team
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

global M_ options_ estim_params_


n = estim_params_.np + ...
    estim_params_.nvn+ ...
    estim_params_.ncx+ ...
    estim_params_.ncn+ ...
    estim_params_.nvx;
nblck = options_.mh_nblck;

MhDirectoryName = CheckPath('metropolis');
load([ MhDirectoryName '/'  M_.fname '_mh_history.mat'])

FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine   = record.KeepedDraws.FirstLine;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));

posterior_kernel_at_the_mode = -Inf;
posterior_mean = zeros(n,1);
posterior_mode = NaN(n,1);
posterior_covariance = zeros(n,n);
offset = 0;

for b=1:nblck
    first_line = FirstLine;
    for n = FirstMhFile:TotalNumberOfMhFiles
        %for b = 1:nblck
        load([ MhDirectoryName '/' M_.fname '_mh' int2str(n) '_blck' int2str(b) '.mat'],'x2','logpo2'); 
        [tmp,idx] = max(logpo2);
        if tmp>posterior_kernel_at_the_mode
            posterior_kernel_at_the_mode = tmp;
            posterior_mode = x2(idx,:);
        end
        [posterior_mean,posterior_covariance,offset] = recursive_moments(posterior_mean,posterior_covariance,x2(first_line:end,:),offset);
        first_line = 1;
    end
end

xparam1 = posterior_mode';
hh = inv(posterior_covariance);
fval = posterior_kernel_at_the_mode;

save([M_.fname '_mh_mode.mat'],'xparam1','hh','fval');