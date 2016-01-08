function bvar_irf(nlags,identification)
% builds IRFs for a bvar model
%
% INPUTS
%    nlags            [integer]     number of lags for the bvar
%    identification   [string]      identification scheme ('Cholesky' or 'SquareRoot')
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2007-2009 Dynare Team
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

if nargin==1
    identification = 'Cholesky';
end

options_ = set_default_option(options_, 'bvar_replic', 2000);

[ny, nx, posterior, prior] = bvar_toolbox(nlags);

S_inv_upper_chol = chol(inv(posterior.S));

% Option 'lower' of chol() not available in old versions of
% Matlab, so using transpose
XXi_lower_chol = chol(posterior.XXi)';

k = ny*nlags+nx;

% Declaration of the companion matrix
Companion_matrix = diag(ones(ny*(nlags-1),1),-ny);

% Number of explosive VAR models sampled
p = 0;

% Initialize a four dimensional array.
sampled_irfs = NaN(ny, ny, options_.irf, options_.bvar_replic);

for draw=1:options_.bvar_replic
    
    % Get a covariance matrix from an inverted Wishart distribution.
    Sigma = rand_inverse_wishart(ny, posterior.df, S_inv_upper_chol);
    Sigma_upper_chol = chol(Sigma);
    Sigma_lower_chol = transpose(Sigma_upper_chol);

    % Get the Autoregressive matrices from a matrix variate distribution.
    Phi = rand_matrix_normal(k, ny, posterior.PhiHat, Sigma_lower_chol, XXi_lower_chol);
    
    % Form the companion matrix.
    Companion_matrix(1:ny,:) = transpose(Phi(1:ny*nlags,:)); 
    
    % All the eigenvalues of the companion matrix have to be on or
    % inside the unit circle to rule out explosive time series.
    test = (abs(eig(Companion_matrix)));
    if any(test>1.0000000000001)
        p = p+1;
    end

    if strcmpi(identification,'Cholesky')
        StructuralMat = Sigma_lower_chol;
    elseif strcmpi(identification,'SquareRoot')
        StructuralMat = sqrtm(Sigma);
    end
    
    % Build the IRFs...
    lags_data = zeros(ny,ny*nlags) ;
    sampled_irfs(:,:,1,draw) = Sigma_lower_chol ;
    lags_data(:,1:ny) = Sigma_lower_chol ;
    for t=2:options_.irf
        sampled_irfs(:,:,t,draw) = lags_data(:,:)*Phi(1:ny*nlags,:) ;
        lags_data(:,ny+1:end) = lags_data(:,1:end-ny) ;
        lags_data(:,1:ny) = sampled_irfs(:,:,t,draw) ;
    end
    
end

if p > 0
    disp(' ')
    disp(['Some of the VAR models sampled from the posterior distribution'])
    disp(['were found to be explosive (' int2str(p) ' samples).'])
    disp(' ')
end

posterior_mean_irfs = mean(sampled_irfs,4);
posterior_variance_irfs = var(sampled_irfs, 1, 4);

sorted_irfs = sort(sampled_irfs,4);
sort_idx = round((0.5 + [-options_.conf_sig, options_.conf_sig, .0]/2) * options_.bvar_replic);

posterior_down_conf_irfs = sorted_irfs(:,:,:,sort_idx(1));
posterior_up_conf_irfs = sorted_irfs(:,:,:,sort_idx(2));
posterior_median_irfs = sorted_irfs(:,:,:,sort_idx(3));   

number_of_columns = fix(sqrt(ny));
number_of_rows = ceil(ny / number_of_columns) ;

% Plots of the IRFs
for shock=1:ny
    figure('Name',['Posterior BVAR Impulse Responses (shock in equation ' int2str(shock) ').']);
    for variable=1:ny
        subplot(number_of_rows,number_of_columns,variable);
        h1 = area(1:options_.irf,squeeze(posterior_up_conf_irfs(shock,variable,:)));
        set(h1,'BaseValue',min([min(posterior_up_conf_irfs(shock,variable,:)),min(posterior_down_conf_irfs(shock,variable,:))]))
        set(h1,'FaceColor',[.9 .9 .9])
        hold on
        h2 = area(1:options_.irf,squeeze(posterior_down_conf_irfs(shock,variable,:)));
        set(h2,'BaseValue',min([min(posterior_up_conf_irfs(shock,variable,:)),min(posterior_down_conf_irfs(shock,variable,:))]))
        set(h2,'FaceColor',[1 1 1])
        plot(1:options_.irf,squeeze(posterior_median_irfs(shock,variable,:)),'-k','linewidth',2)
        axis tight
        hold off
    end
end

% Save intermediate results
DirectoryName = [ M_.fname '/bvar_irf' ];
if ~isdir(DirectoryName)
    mkdir('.',DirectoryName);
end
save([ DirectoryName '/simulations.mat'], 'sampled_irfs');

% Save results in oo_
for i=1:ny
    shock_name = options_.varobs(i, :);
    for j=1:ny
        variable_name = options_.varobs(j, :);
        eval(['oo_.bvar.irf.Mean.' variable_name '.' shock_name ' = posterior_mean_irfs(' int2str(j) ',' int2str(i) ',:);'])
        eval(['oo_.bvar.irf.Median.' variable_name '.' shock_name ' = posterior_median_irfs(' int2str(j) ',' int2str(i) ',:);'])
        eval(['oo_.bvar.irf.Var.' variable_name '.' shock_name ' = posterior_variance_irfs(' int2str(j) ',' int2str(i) ',:);'])
        eval(['oo_.bvar.irf.Upper_bound.' variable_name '.' shock_name ' = posterior_up_conf_irfs(' int2str(j) ',' int2str(i) ',:);'])
        eval(['oo_.bvar.irf.Lower_bound.' variable_name '.' shock_name ' = posterior_down_conf_irfs(' int2str(j) ',' int2str(i) ',:);'])
    end
end