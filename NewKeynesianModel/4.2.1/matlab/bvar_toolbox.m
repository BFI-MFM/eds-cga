function [ny, nx, posterior, prior, forecast_data] = bvar_toolbox(nlags)
%function [ny, nx, posterior, prior, forecast_data] = bvar_toolbox(nlags)
% bvar_toolbox  Routines shared between BVAR methods
% Computes several things for the estimations of a BVAR(nlags)
%
% INPUTS:
%     nlags: number of lags
%
% OUTPUTS:
%    ny:            number of endogenous variables
%    nx:            number of exogenous variables (equal to zero, or one if a
%                   constant term is included)
%    posterior:     a structure describing the posterior distribution (which is
%                   normal-Inverse-Wishart)
%                   Its fields are:
%                   - df: degrees of freedom of the inverse-Wishart distribution
%                   - S: matrix parameter for the inverse-Wishart distribution
%                   - XXi: first component of the VCV of the matrix-normal
%                     distribution (the other one being drawn from the
%                     inverse-Wishart)
%                   - PhiHat: mean of the matrix-normal distribution
%    prior:         a structure describing the prior distribution
%                   Its fields are the same than for the posterior
%    forecast_data: a structure containing data useful for forecasting
%                   Its fields are:
%                   - initval: a nlags*ny matrix containing the "nlags" last
%                     observations of the sample (i.e. before options_.nobs)
%                   - xdata: a matrix containing the future exogenous for
%                     forecasting, of size options_.forecast*nx (actually only
%                     contains "1" values for the constant term if nx ~= 0)
%                   - realized_val: only non-empty if options_.nobs doesn't point
%                     to the end of sample    
%                     In that case, contains values of endogenous variables after
%                     options_.nobs and up to the end of the sample
%                   - realized_xdata: contains values of exogenous variables after
%                     options_.nobs and up to the end of the sample (actually only
%                     contains "1" values for the constant term if nx ~= 0)
%
% SPECIAL REQUIREMENTS:
%    This function uses the following Dynare options:
%    - datafile, first_obs, varobs, xls_sheet, xls_range, nobs, presample
%    - bvar_prior_{tau,decay,lambda,mu,omega,flat,train}

% Copyright (C) 2003-2007 Christopher Sims
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

global options_

% Load dataset
dataset = read_variables(options_.datafile, options_.varobs, [], options_.xls_sheet, options_.xls_range);
options_ = set_default_option(options_, 'nobs', size(dataset,1)-options_.first_obs+1);

% Parameters for prior
options_ = set_default_option(options_, 'bvar_prior_tau', 3);
options_ = set_default_option(options_, 'bvar_prior_decay', 0.5);
options_ = set_default_option(options_, 'bvar_prior_lambda', 5);
options_ = set_default_option(options_, 'bvar_prior_mu', 2);
options_ = set_default_option(options_, 'bvar_prior_omega', 1);
options_ = set_default_option(options_, 'bvar_prior_flat', 0);
options_ = set_default_option(options_, 'bvar_prior_train', 0);

if options_.first_obs + options_.presample <= nlags
    error('first_obs+presample should be > nlags (for initializing the VAR)')
end

train = options_.bvar_prior_train;

if options_.first_obs + options_.presample - train <= nlags
    error('first_obs+presample-train should be > nlags (for initializating the VAR)')
end

idx = options_.first_obs+options_.presample-train-nlags:options_.first_obs+options_.nobs-1;

% Prepare dataset
if options_.loglinear && ~options_.logdata
    dataset = log(dataset);
end
if options_.prefilter
    dataset(idx,:) = dataset(idx,:) - ones(length(idx),1)*mean(dataset(idx,:));
end

mnprior.tight = options_.bvar_prior_tau;
mnprior.decay = options_.bvar_prior_decay;

% Use only initializations lags for the variance prior
vprior.sig = std(dataset(options_.first_obs+options_.presample-nlags:options_.first_obs+options_.presample-1,:))';
vprior.w = options_.bvar_prior_omega;

lambda = options_.bvar_prior_lambda;
mu = options_.bvar_prior_mu;
flat = options_.bvar_prior_flat;

ny = size(dataset, 2);
if options_.prefilter || options_.noconstant
    nx = 0;
else
    nx = 1;
end

[ydum, xdum, pbreaks] = varprior(ny, nx, nlags, mnprior, vprior);

ydata = dataset(idx, :);
T = size(ydata, 1);
xdata = ones(T,nx);

% Posterior density
var = rfvar3([ydata; ydum], nlags, [xdata; xdum], [T; T+pbreaks], lambda, mu);
Tu = size(var.u, 1);

posterior.df = Tu - ny*nlags - nx - flat*(ny+1);
posterior.S = var.u' * var.u;
posterior.XXi = var.xxi;
posterior.PhiHat = var.B;

% Prior density
Tp = train + nlags;
if nx
    xdata = xdata(1:Tp, :);
else
    xdata = [];
end
varp = rfvar3([ydata(1:Tp, :); ydum], nlags, [xdata; xdum], [Tp; Tp + pbreaks], lambda, mu);
Tup = size(varp.u, 1);

prior.df = Tup - ny*nlags - nx - flat*(ny+1);
prior.S = varp.u' * varp.u;
prior.XXi = varp.xxi;
prior.PhiHat = varp.B;

if prior.df < ny
    error('Too few degrees of freedom in the inverse-Wishart part of prior distribution. You should increase training sample size.')
end

% Add forecast informations
if nargout >= 5
    forecast_data.xdata = ones(options_.forecast, nx);
    forecast_data.initval = ydata(end-nlags+1:end, :);
    if options_.first_obs + options_.nobs <= size(dataset, 1)
        forecast_data.realized_val = dataset(options_.first_obs+options_.nobs:end, :);
        forecast_data.realized_xdata = ones(size(forecast_data.realized_val, 1), nx);
    else
        forecast_data.realized_val = [];
    end
end


function [ydum,xdum,breaks]=varprior(nv,nx,lags,mnprior,vprior)
%function [ydum,xdum,breaks]=varprior(nv,nx,lags,mnprior,vprior)
% ydum, xdum:   dummy observation data that implement the prior
% breaks:       vector of points in the dummy data after which new dummy obs's start
%                   Set breaks=T+[0;breaks], ydata=[ydata;ydum], xdum=[xdata;xdum], where 
%                   actual data matrix has T rows, in preparing input for rfvar3
% nv,nx,lags: VAR dimensions
% mnprior.tight:Overall tightness of Minnesota prior
% mnprior.decay:Standard deviations of lags shrink as lag^(-decay)
% vprior.sig:   Vector of prior modes for diagonal elements of r.f. covariance matrix
% vprior.w:     Weight on prior on vcv.  1 corresponds to "one dummy observation" weight
%                   Should be an integer, and will be rounded if not.  vprior.sig is needed
%                   to scale the Minnesota prior, even if the prior on sigma is not used itself.
%                   Set vprior.w=0 to achieve this.
% Note:         The original Minnesota prior treats own lags asymmetrically, and therefore
%                   cannot be implemented entirely with dummy observations.  It is also usually
%                   taken to include the sum-of-coefficients and co-persistence components
%                   that are implemented directly in rfvar3.m.  The diagonal prior on v, combined
%                   with sum-of-coefficients and co-persistence components and with the unit own-first-lag
%                   prior mean generates larger prior variances for own than for cross-effects even in 
%                   this formulation, but here there is no way to shrink toward a set of unconstrained 
%                   univariate AR's.

% Original file downloaded from:
% http://sims.princeton.edu/yftp/VARtools/matlab/varprior.m

if ~isempty(mnprior)
    xdum = zeros(lags+1,nx,lags,nv);
    ydum = zeros(lags+1,nv,lags,nv);
    for il = 1:lags
        ydum(il+1,:,il,:) = il^mnprior.decay*diag(vprior.sig);
    end
    ydum(1,:,1,:) = diag(vprior.sig);
    ydum = mnprior.tight*reshape(ydum,[lags+1,nv,lags*nv]);
    ydum = flipdim(ydum,1);
    xdum = mnprior.tight*reshape(xdum,[lags+1,nx,lags*nv]);
    xdum = flipdim(xdum,1);
    breaks = (lags+1)*[1:(nv*lags)]';
    lbreak = breaks(end);
else
    ydum = [];
    xdum = [];
    breaks = [];
    lbreak = 0;
end
if ~isempty(vprior) && vprior.w>0
    ydum2 = zeros(lags+1,nv,nv);
    xdum2 = zeros(lags+1,nx,nv);
    ydum2(end,:,:) = diag(vprior.sig);
    for i = 1:vprior.w
        ydum = cat(3,ydum,ydum2);
        xdum = cat(3,xdum,xdum2);
        breaks = [breaks;(lags+1)*[1:nv]'+lbreak];
        lbreak = breaks(end);
    end
end
dimy = size(ydum);
ydum = reshape(permute(ydum,[1 3 2]),dimy(1)*dimy(3),nv);
xdum = reshape(permute(xdum,[1 3 2]),dimy(1)*dimy(3),nx);
breaks = breaks(1:(end-1));


function var=rfvar3(ydata,lags,xdata,breaks,lambda,mu)
%function var=rfvar3(ydata,lags,xdata,breaks,lambda,mu)
% This algorithm goes for accuracy without worrying about memory requirements.
% ydata:   dependent variable data matrix
% xdata:   exogenous variable data matrix
% lags:    number of lags
% breaks:  rows in ydata and xdata after which there is a break.  This allows for
%          discontinuities in the data (e.g. war years) and for the possibility of
%          adding dummy observations to implement a prior.  This must be a column vector.
%          Note that a single dummy observation becomes lags+1 rows of the data matrix,
%          with a break separating it from the rest of the data.  The function treats the 
%          first lags observations at the top and after each "break" in ydata and xdata as
%          initial conditions. 
% lambda:  weight on "co-persistence" prior dummy observations.  This expresses
%          belief that when data on *all* y's are stable at their initial levels, they will
%          tend to persist at that level.  lambda=5 is a reasonable first try.  With lambda<0,
%          constant term is not included in the dummy observation, so that stationary models
%          with means equal to initial ybar do not fit the prior mean.  With lambda>0, the prior
%          implies that large constants are unlikely if unit roots are present.
% mu:      weight on "own persistence" prior dummy observation.  Expresses belief
%          that when y_i has been stable at its initial level, it will tend to persist
%          at that level, regardless of the values of other variables.  There is
%          one of these for each variable.  A reasonable first guess is mu=2.
%      The program assumes that the first lags rows of ydata and xdata are real data, not dummies.
%      Dummy observations should go at the end, if any.  If pre-sample x's are not available,
%      repeating the initial xdata(lags+1,:) row or copying xdata(lags+1:2*lags,:) into 
%      xdata(1:lags,:) are reasonable subsititutes.  These values are used in forming the
%      persistence priors.

% Original file downloaded from:
% http://sims.princeton.edu/yftp/VARtools/matlab/rfvar3.m

[T,nvar] = size(ydata);
nox = isempty(xdata);
if ~nox
    [T2,nx] = size(xdata);
else
    T2 = T;
    nx = 0;
    xdata = zeros(T2,0);
end
% note that x must be same length as y, even though first part of x will not be used.
% This is so that the lags parameter can be changed without reshaping the xdata matrix.
if T2 ~= T, error('Mismatch of x and y data lengths'),end
if nargin < 4
    nbreaks = 0;
    breaks = [];
else
    nbreaks = length(breaks);
end
breaks = [0;breaks;T];
smpl = [];
for nb = 1:nbreaks+1
    smpl = [smpl;[breaks(nb)+lags+1:breaks(nb+1)]'];
end
Tsmpl = size(smpl,1);
X = zeros(Tsmpl,nvar,lags);
for is = 1:length(smpl)
    X(is,:,:) = ydata(smpl(is)-(1:lags),:)';
end
X = [X(:,:) xdata(smpl,:)];
y = ydata(smpl,:);
% Everything now set up with input data for y=Xb+e 

% Add persistence dummies
if lambda ~= 0 || mu > 0
    ybar = mean(ydata(1:lags,:),1);
    if ~nox
        xbar = mean(xdata(1:lags,:),1);
    else
        xbar = [];
    end
    if lambda ~= 0
        if lambda>0
            xdum = lambda*[repmat(ybar,1,lags) xbar];
        else
            lambda = -lambda;
            xdum = lambda*[repmat(ybar,1,lags) zeros(size(xbar))];
        end
        ydum = zeros(1,nvar);
        ydum(1,:) = lambda*ybar;
        y = [y;ydum];
        X = [X;xdum];
    end
    if mu>0
        xdum = [repmat(diag(ybar),1,lags) zeros(nvar,nx)]*mu;
        ydum = mu*diag(ybar);
        X = [X;xdum];
        y = [y;ydum];
    end
end

% Compute OLS regression and residuals
[vl,d,vr] = svd(X,0);
di = 1./diag(d);
B = (vr.*repmat(di',nvar*lags+nx,1))*vl'*y;
u = y-X*B;
xxi = vr.*repmat(di',nvar*lags+nx,1);
xxi = xxi*xxi';

var.B = B;
var.u = u;
var.xxi = xxi;
