%
% SYNOPSIS
% 
% r = dynare_simul(name, shocks)
% r = dynare_simul(name, prefix, shocks)
% r = dynare_simul(name, shocks, start)
% r = dynare_simul(name, prefix, shocks, start)
%
%     name     name of MAT-file produced by dynare++
%     prefix   prefix of variables in the MAT-file
%     shocks   matrix of shocks
%     start    zero period value
%
% Note that this file requires the dynare_simul_ DLL to be in the path.
% This DLL is distributed with Dynare, under the mex/matlab or mex/octave
% subdirectory.
%
% SEMANTICS
%
% The command reads a decision rule from the MAT-file having the given
% prefix. Then it starts simulating the decision rule with zero time value
% equal to the given start. It uses the given shocks for the simulation. If
% the start is not given, the state about which the decision rule is
% centralized is taken (called fix point, or stochastic steady state, take
% your pick).
%
%     prefix   Use the prefix with which you called dynare++, the default
%              prefix in dynare++ is 'dyn'.
%     shocks   Number of rows must be a number of exogenous shocks,
%              number of columns gives the number of simulated
%              periods. NaNs and Infs in the matrix are substitued by
%              draws from the normal distribution using the covariance
%              matrix given in the model file.
%     start    Vector of endogenous variables in the ordering given by
%              <prefix>_vars.
%
% Seed for random generator is derived from calling rand(1,1). Therefore,
% seeding can be controlled with rand('state') and rand('state',some_seed).
%
% EXAMPLES
%
% All examples suppose that the prefix is 'dyn' and that your_model.mat
% has been loaded into Matlab.
%
% 1. response to permanent negative shock to the third exo var EPS3 for
%    100 periods
%
%       shocks = zeros(4,100); % 4 exogenous variables in the model
%       shocks(dyn_i_EPS3,:) = -0.1; % the permanent shock to EPS3
%       r = dynare_simul('your_model.mat',shocks);
%
% 2. one stochastic simulation for 100 periods
%
%       shocks = zeros(4,100)./0; % put NaNs everywhere
%       r = dynare_simul('your_model.mat',shocks);
%
% 3. one stochastic simulation starting at 75% undercapitalized economy
%
%       shocks = zeros(4,100)./0; % put NaNs everywhere
%       ystart = dyn_ss; % get copy of DR fix point
%       ystart(dyn_i_K) = 0.75*dyn_ss(dyn_i_K); % scale down the capital
%       r = dynare_simul('your_model.mat',shocks,ystart);
%
% 
% SEE ALSO
%
%   "DSGE Models with Dynare++. A Tutorial.", Ondra Kamenik, 2005

% Copyright (C) 2005-2011, Ondra Kamenik

function r = dynare_simul(varargin)

if exist('dynare_simul_') ~= 3
    error('Can''t find dynare_simul_ DLL in the path. The simplest way to add it is to run Dynare once in this session.')
end

% get the file name and load data
fname = varargin{1};
eval(['load ' fname]);

% set prefix, shocks, ystart
if ischar(varargin{2})
  prefix = varargin{2};
  if length(varargin) == 3
    shocks = varargin{3};
    ystart = NaN;
  elseif length(varargin) == 4
    shocks = varargin{3};
    ystart = varargin{4};
  else
    error('Wrong number of parameters.');
  end
else
  prefix = 'dyn';
  if length(varargin) == 2
    shocks = varargin{2};
    ystart = NaN;
  elseif length(varargin) == 3
    shocks = varargin{2};
    ystart = varargin{3};
  else
    error('Wrong number of parameters.');
  end
end

% load all needed variables but prefix_g_*
if (exist([prefix '_nstat']))
  nstat = eval([prefix '_nstat']);
else
  error(['Could not find variable ' prefix '_nstat in workspace']);
end
if (exist([prefix '_npred']))
  npred = eval([prefix '_npred']);
else
  error(['Could not find variable ' prefix '_npred in workspace']);
end
if (exist([prefix '_nboth']))
  nboth = eval([prefix '_nboth']);
else
  error(['Could not find variable ' prefix '_nboth in workspace']);
end
if (exist([prefix '_nforw']))
  nforw = eval([prefix '_nforw']);
else
  error(['Could not find variable ' prefix '_nforw in workspace']);
end
if (exist([prefix '_ss']))
  ss = eval([prefix '_ss']);
else
  error(['Could not find variable ' prefix '_ss in workspace']);
end
if (exist([prefix '_vcov_exo']))
  vcov_exo = eval([prefix '_vcov_exo']);
else
  error(['Could not find variable ' prefix '_vcov_exo in workspace']);
end
nexog = size(vcov_exo,1);

if isnan(ystart)
  ystart = ss;
end

% newer version of dynare++ doesn't return prefix_g_0, we make it here if
% it does not exist in workspace
g_zero = [prefix '_g_0'];
if (~ exist(g_zero))
  eval([ g_zero '= zeros(nstat+npred+nboth+nforw,1);']);
end

% make derstr a string of comma seperated existing prefix_g_*
derstr = [',' g_zero];
order = 1;
cont = 1;
while cont == 1
  g_ord = [prefix '_g_' num2str(order)];
  if (exist(g_ord))
    derstr = [derstr ',' g_ord];
    order = order + 1;
  else
    cont = 0;
  end
end

% set seed
seed = ceil(10000*rand(1,1));

% call dynare_simul_
command = ['[err,r]=dynare_simul_(' num2str(order-1) ',nstat,npred,nboth,nforw,' ...
           'nexog,ystart,shocks,vcov_exo,seed,ss' derstr ');'];
eval(command);

if err
    error('Simulation failed')
end
