function bvar_density(maxnlags)
% function bvar_density(maxnlags)
% computes the density of a bayesian var
%
% INPUTS
%    maxnlags:      maximum number of lags in the bvar
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2007 Christopher Sims
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

for nlags = 1:maxnlags
    [ny, nx, posterior, prior] = bvar_toolbox(nlags);
    
    posterior_int = matrictint(posterior.S, posterior.df, posterior.XXi);
    prior_int = matrictint(prior.S, prior.df, prior.XXi);
    
    lik_nobs = posterior.df - prior.df;
    
    log_dnsty = posterior_int - prior_int - 0.5*ny*lik_nobs*log(2*pi);
    
    disp(' ')
    fprintf('The marginal log density of the BVAR(%g) model is equal to %10.4f\n', ...
            nlags, log_dnsty);
    disp(' ')
end


function w = matrictint(S, df, XXi)
% Computes the log of the integral of the kernel of the PDF of a
% normal-inverse-Wishart distribution.
%
% S:   parameter of inverse-Wishart distribution
% df:  number of degrees of freedom of inverse-Wishart distribution
% XXi: first component of VCV matrix of matrix-normal distribution
% 
% Computes the integral over (Phi, Sigma) of:
%
% det(Sigma)^(-k/2)*exp(-0.5*Tr((Phi-PhiHat)'*(XXi)^(-1)*(Phi-PhiHat)*Sigma^(-1)))*
% det(Sigma)^((df+ny+1)/2)*exp(-0.5*Tr(Sigma^(-1)*S))
%
% (where k is the dimension of XXi and ny is the dimension of S and
% Sigma)

% Original file downloaded from:
% http://sims.princeton.edu/yftp/VARtools/matlab/matrictint.m

k=size(XXi,1);
ny=size(S,1);
[cx,p]=chol(XXi);
[cs,q]=chol(S);

if any(diag(cx)<100*eps)
    error('singular XXi')
end
if any(diag(cs<100*eps))
    error('singular S')
end

% Matrix-normal component
w1 = 0.5*k*ny*log(2*pi)+ny*sum(log(diag(cx)));

% Inverse-Wishart component
w2 = -df*sum(log(diag(cs))) + 0.5*df*ny*log(2) + ny*(ny-1)*0.25*log(pi) + ggammaln(ny, df);

w = w1 + w2;

function lgg = ggammaln(m, df)
if df <= (m-1)
    error('too few df in ggammaln')
else
    garg = 0.5*(df+(0:-1:1-m));
    lgg = sum(gammaln(garg));
end
